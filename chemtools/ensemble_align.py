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
parser.add_argument('--smarts', help='Ignore MCS determination, use an user-provided SMARTS string')

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

if args.smarts:
    log.write(f'User-provided SMARTS string {args.smarts}\n')
    pattern = Chem.MolFromSmarts(args.smarts)

else:
    log.write('Determining MCS (Maximum Common Substructure)\n')
    mcs = rdFMCS.FindMCS(all_mols)
    pattern = Chem.MolFromSmarts(mcs.smartsString)
    log.write('MCS was determined using all atoms\n')
    log.write(f"MCS has {mcs.numAtoms} atoms\n")
    log.write(f"MCS SMARTS Pattern is {mcs.smartsString}\n\n")
    
log.write('Removing hydrogens from the SMARTS patterns\n')
pattern = Chem.RemoveAllHs(pattern)
mcs_heavy = Chem.MolToSmarts(pattern)
log.write(f'MCS SMARTS Pattern (heavy-atoms) is {mcs_heavy}\n\n')

ref_match = ref_mol.GetSubstructMatch(pattern)

log.write(f'Aligning ensemble {args.ensemble} to reference {args.reference}\n')
calculated_rmsd = {}
for i, probe_mol in enumerate(ens_mols):
    mv = probe_mol.GetSubstructMatch(pattern)
    rmsd = rdMolAlign.AlignMol(probe_mol, ref_mol, atomMap=list(zip(mv,ref_match)))
    calculated_rmsd[f'{ens_header[i]}'] = rmsd

log.write('Writing aligned ensemble to file\n')
with open(f"{args.ensemble.rsplit('.')[0]}_aligned.xyz", mode='w') as f:
    for mol in ens_mols:
        f.write(Chem.MolToXYZBlock(mol))

calculated_rmsd = pd.Series(calculated_rmsd, name='rmsd', dtype=float)
if args.reset_names:
    basename = args.ensemble
    basename = basename.rsplit('.',1)[0]
    calculated_rmsd.index = [f"{basename}_{i+1}.xyz" for i in range(len(calculated_rmsd))] # type: ignore
if args.sort:
    calculated_rmsd.sort_values(inplace=True)
calculated_rmsd.index.name='file'

log.write('Writing calculated RMSD to file\n')
calculated_rmsd.to_csv('rmsd.csv', float_format='%.5f')

log.write('Finished\n\n')
log.close()