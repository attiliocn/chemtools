#! /usr/bin/env python3

import argparse
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd
from modules import rdkitutils, xyzutils
from rdkit import Chem
from rdkit.Chem import rdMolAlign, rdMolTransforms, rdDetermineBonds, rdFMCS

parser = argparse.ArgumentParser()
parser.add_argument('reference', help='Reference Molecule (XYZ format)')
parser.add_argument('ensemble', help='Ensemble of molecules (XYZ format)')
parser.add_argument('--mcs', action='store_true', help='Try to determine the Maximum Common Substructure (MCS) and align molecules using the resulting MCS')
parser.add_argument('--rmsd-all', action='store_true', help='Use all-atom RMSD calculation. Default is False')
args = parser.parse_args()

# read the reference molecule
ref_elements, ref_coordinates = xyzutils.read_xyz_file(args.reference)
ref_mol = rdkitutils.convert_coordinates_to_mols(ref_elements, ref_coordinates)

# read the ensemble to be aligned
ensemble = xyzutils.read_xyz_ensemble(args.ensemble)
ens_elements = [_['elements'] for _ in ensemble.values()]
ens_coordinates = [_['coordinates'] for _ in ensemble.values()]
ens_header = [_['header'] for _ in ensemble.values()]
ens_mols = [rdkitutils.convert_coordinates_to_mols(ele, coords) for ele, coords in zip(ens_elements, ens_coordinates)]
numMols = len(ens_mols)

# combine the reference molecule with the ensemble molecule
full_ensemble = [ref_mol] + ens_mols

if not args.rmsd_all:
    # remove Hs from all molecules in the full ensemble
    full_ensemble = [Chem.RemoveAllHs(mol) for mol in full_ensemble]
    for mol in full_ensemble:
        rdDetermineBonds.DetermineConnectivity(mol)

if args.mcs:
    mcs = rdFMCS.FindMCS(full_ensemble)
    pattern = Chem.MolFromSmarts(mcs.smartsString)
    print('MCS was determined using all atoms\n')
    print(f"MCS has {mcs.numAtoms} atoms\n")
    print(f"MCS SMARTS Pattern is {mcs.smartsString}\n\n")
    full_ensemble[0] = Chem.MolFromSmarts(mcs.smartsString)

# get the transformation matrices
def get_best_align_transf(args):
    reference, probe = args

    atomMap = []
    # try to map PROBE atoms into REFERENCE atoms
    numAtoms = probe.GetNumAtoms()
    substructMatches = probe.GetSubstructMatches(reference, uniquify=False)
    for refAtom in substructMatches:
        atomMap.append(list(zip(range(numAtoms), refAtom)))
    
    # if fails, try to map REFERENCE into PROBE
    if not substructMatches:
        numAtoms = reference.GetNumAtoms()
        substructMatches = reference.GetSubstructMatches(probe, uniquify=False)
        for probeAtom in substructMatches:
            atomMap.append(list(zip(probeAtom,range(numAtoms))))
     
    rmsd, R, atomMap = rdMolAlign.GetBestAlignmentTransform(reference,probe,map=atomMap)
    return rmsd, R, atomMap
tasks = [(full_ensemble[0], probe) for probe in full_ensemble[1:]]
with Pool(processes=os.cpu_count()) as pool:
    results = pool.map(get_best_align_transf, tasks)

# parse and export results
rmsd_all = [result[0] for result in results]
transform_matrices = [result[1] for result in results]

for i in range(numMols):
    rdMolTransforms.TransformConformer(ens_mols[i].GetConformer(), transform_matrices[i])

with open(f"{args.ensemble.rsplit('.')[0]}_aligned.xyz", mode='w') as f:
    for mol in ens_mols:
        f.write(Chem.MolToXYZBlock(mol, confId=0))

calculated_rmsd = {k:v for k,v in zip(ens_header,rmsd_all)}
calculated_rmsd = pd.Series(calculated_rmsd, name='rmsd', dtype=float)
calculated_rmsd.to_csv('rmsd.csv', float_format='%.4f')