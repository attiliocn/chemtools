#! /usr/bin/env python3

import argparse
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import rdFMCS, rdMolAlign

from modules.xyzutils import read_xyz_ensemble, read_xyz_file
from modules.rdkitutils import convert_coordinates_to_mols

parser = argparse.ArgumentParser()
parser.add_argument('reference', help='Reference Molecule (XYZ format)')
parser.add_argument('ensemble', help='Ensemble of molecules (XYZ format)')
parser.add_argument('--sort', action='store_true', help='Sort the log file by increasing RMSD to the reference')
parser.add_argument('--reset-names', action='store_true', help='Rename the ensemble structures')
args = parser.parse_args()

log = open('align.log', mode='w', buffering=1)

ref_elements, ref_coordinates = read_xyz_file(args.reference)
ref_mol = convert_coordinates_to_mols(ref_elements, ref_coordinates)

ensemble = read_xyz_ensemble(args.ensemble)
ens_elements = [_['elements'] for _ in ensemble.values()]
ens_coordinates = [_['coordinates'] for _ in ensemble.values()]
ens_header = [_['header'] for _ in ensemble.values()]
ens_mols = [convert_coordinates_to_mols(ele, coords) for ele, coords in zip(ens_elements, ens_coordinates)]

all_mols = np.concatenate([[ref_mol],ens_mols])
mcs = rdFMCS.FindMCS(all_mols)
pattern = Chem.MolFromSmarts(mcs.smartsString)

pattern = Chem.RemoveAllHs(pattern)
mcs_heavy = Chem.MolToSmarts(pattern)

ref_match = ref_mol.GetSubstructMatch(pattern)

calculated_rmsd = {}
for i, probe_mol in enumerate(ens_mols):
    mv = probe_mol.GetSubstructMatch(pattern)
    rmsd = rdMolAlign.AlignMol(probe_mol, ref_mol, atomMap=list(zip(mv,ref_match)))
    calculated_rmsd[f'{ens_header[i]}'] = rmsd

with open(f"{args.ensemble.rsplit('.')[0]}_aligned.xyz", mode='w') as f:
    for mol in ens_mols:
        f.write(Chem.MolToXYZBlock(mol))

calculated_rmsd = pd.Series(calculated_rmsd)
if args.reset_names:
    basename = args.ensemble
    basename = basename.rsplit('.',1)[0]
    calculated_rmsd.index = [f"{basename}_{i+1}.xyz" for i in range(len(calculated_rmsd))]
if args.sort:
    calculated_rmsd.sort_values(inplace=True)

for k,v in calculated_rmsd.items():
    log.write(f"{k},{v:.5f}\n")

log.close()