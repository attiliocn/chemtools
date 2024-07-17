#!/usr/bin/env python3

import argparse
from modules import xyzutils, rdkitutils
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign
import os

parser = argparse.ArgumentParser()
parser.add_argument('ensemble', help='XYZ Ensemble File')
parser.add_argument('--max-matches', '-c', type=int, default=1000000, help='Maximum number of symmetric substructures to consider. Default is 1e6')
args = parser.parse_args()

basename, extension = args.ensemble.rsplit('.', 1)

ensemble = xyzutils.read_xyz_ensemble(args.ensemble)
numconfs = len(ensemble)
ens_elements = [_['elements'] for _ in ensemble.values()]
ens_coordinates = [_['coordinates'] for _ in ensemble.values()]
ens_header = [_['header'] for _ in ensemble.values()]
ens_mols = [rdkitutils.convert_coordinates_to_mols(ele, coords, removeHs=True) for ele, coords in zip(ens_elements, ens_coordinates)]

mol_with_conformers = ens_mols[0]
for mol in ens_mols[1:]:
    conformer = mol.GetConformer()
    mol_with_conformers.AddConformer(conformer,assignId=True)

atomMap = rdkitutils.get_maximum_substructure_matches(ens_mols, max_matches=args.max_matches)
conformers_rmsd = rdMolAlign.GetAllConformerBestRMS(
    mol=mol_with_conformers,
    map=atomMap,
    symmetrizeConjugatedTerminalGroups=True,
    numThreads=os.cpu_count(),
)

rmsd_distance_matrix = np.zeros(numconfs**2).reshape(numconfs,numconfs)
k = 0
for i in range(numconfs):
    for j in range(i):
        rmsd_distance_matrix[i,j] = conformers_rmsd[k]
        k+=1
np.savetxt(f"{basename}.csv", rmsd_distance_matrix, delimiter=",")

