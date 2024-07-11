import os
from multiprocessing import Pool

import numpy as np
from modules.xyzutils import build_xyz_file
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, rdMolAlign

def convert_coordinates_to_mols(elements, coordinates):
    xyz = build_xyz_file(elements, coordinates)
    mol = Chem.rdmolfiles.MolFromXYZBlock(xyz)
    rdDetermineBonds.DetermineConnectivity(mol)
    return mol

def rmsd(probe_mol,ref_mol):
    rmsd = rdMolAlign.GetBestRMS(probe_mol, ref_mol)
    return rmsd
def calculate_rmsd(args):
    i, j, mols = args
    return i, j, rmsd(mols[i], mols[j])
def rmsd_matrix_parallel(mols:list):
    num_matrices = len(mols)
    rmsd_matrix = np.zeros((num_matrices, num_matrices))
    tasks = [(i, j, mols) for i in range(num_matrices) for j in range(i)]
    with Pool(processes=os.cpu_count()) as pool:
        results = pool.map(calculate_rmsd, tasks)
    for i, j, rmsd_value in results:
        rmsd_matrix[i, j] = rmsd_value
    return rmsd_matrix