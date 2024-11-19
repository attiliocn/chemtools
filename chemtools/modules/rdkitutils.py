import os
from multiprocessing import Pool, shared_memory

import numpy as np
from chemtools.modules import xyzutils
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, rdMolAlign

def convert_coordinates_to_mols(elements, coordinates, removeHs=False):
    xyz = xyzutils.build_xyz_file(elements, coordinates)
    mol = Chem.rdmolfiles.MolFromXYZBlock(xyz)
    if removeHs:
        mol = Chem.RemoveAllHs(mol)
    rdDetermineBonds.DetermineConnectivity(mol)
    return mol

def get_maximum_substructure_matches(mols, max_matches=1000):
    '''
    NOTE: this function
    does not work if molecules within the ensemble
    do not have the same atom numbers
    '''
    matches = []
    for i in range(len(mols)):
        for j in range(i):
            match = mols[i].GetSubstructMatches(mols[j], uniquify=False, maxMatches=max_matches)
            if len(match) > len(matches):
                matches = match
        if len(matches) >= max_matches:
            break
    atomMap = []
    for match in matches:
        atomMap.append(list(zip(range(mols[0].GetNumAtoms()), match)))
    return atomMap

def share_mols_to_ram(mols:list)->str:
    # convert mols from rdkit format to bytes
    mols_bytes = [mol.ToBinary() for mol in mols]
    # share the bytes-converted molecules to RAM
    mols_shared = shared_memory.ShareableList(mols_bytes)
    print(f"Size of the shared array {mols_shared.shm.size/1e6:.2f} MB")
    return mols_shared.shm.name

def get_substructure_matches(args):
    mols_name, i, j,  kwargs = args
    shared_mols = shared_memory.ShareableList(name=mols_name)
    mol_1 = Chem.Mol(shared_mols[i])
    mol_2 = Chem.Mol(shared_mols[j])
    substructures = mol_1.GetSubstructMatches(mol_2, **kwargs)
    shared_mols.shm.close()
    return substructures
def get_symmetric_substructures(mols:list, maxMatches:int=1000):
    # expose the mol list to a shared RAM 
    shared_mols_name = share_mols_to_ram(mols)

    tasks = []
    for i in range(len(mols)):
        for j in range(i):
            tasks.append(
                (shared_mols_name, i, j, {'maxMatches': maxMatches, 'uniquify': False})
            )
    with Pool(processes=os.cpu_count()) as pool:
        matches = []
        results = pool.imap_unordered(get_substructure_matches, tasks)
        for substructures in results:
            if len(substructures) > len(matches):
                matches = substructures
                print(f'Current substructure count is {len(matches)}')
            elif len(matches) >= maxMatches:
                print(f'The total substructures count hit MaxMatches={maxMatches}. The paralell execution will terminate.')
                break
        
        shared_memory.ShareableList(name=shared_mols_name).shm.unlink()

    atomMap = []
    for match in matches:
        atomMap.append(list(zip(range(mols[0].GetNumAtoms()), match)))
    return atomMap

def rmsd(probe_mol,ref_mol,atomMap=None):
    rmsd = rdMolAlign.GetBestRMS(probe_mol, ref_mol, map=atomMap)
    return rmsd

def calculate_rmsd(args):
    i, j, mols, atomMap = args
    return i, j, rmsd(mols[i], mols[j], atomMap=atomMap)
def rmsd_matrix_parallel(mols:list, atomMap=None):
    num_matrices = len(mols)
    rmsd_matrix = np.zeros((num_matrices, num_matrices))
    tasks = [(i, j, mols, atomMap) for i in range(num_matrices) for j in range(i)]
    with Pool(processes=os.cpu_count()) as pool:
        results = pool.map(calculate_rmsd, tasks)
    for i, j, rmsd_value in results:
        rmsd_matrix[i, j] = rmsd_value
    return rmsd_matrix