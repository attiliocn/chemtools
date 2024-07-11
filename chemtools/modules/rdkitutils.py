import os
from multiprocessing import Pool

import numpy as np
from modules.xyzutils import build_xyz_file
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, rdMolAlign

def convert_coordinates_to_mols(elements, coordinates):
    xyz = build_xyz_file(elements, coordinates)
    mol = Chem.rdmolfiles.MolFromXYZBlock(xyz)
    mol = Chem.RemoveAllHs(mol)
    rdDetermineBonds.DetermineConnectivity(mol)
    return mol

def get_maximum_substructure_matches(mols, max_matches=1000):
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

def rmsd(probe_mol,ref_mol, atomMap=None):
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