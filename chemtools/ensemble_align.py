#! /usr/bin/env python3

import argparse

import numpy as np
import pandas as pd
from modules import rdkitutils, xyzutils
from rdkit import Chem

parser = argparse.ArgumentParser()
parser.add_argument('reference', help='Reference Molecule (XYZ format)')
parser.add_argument('ensemble', help='Ensemble of molecules (XYZ format)')
parser.add_argument('--removeH-struc', action='store_true', help='Remove all hydrogens from all molecules. Default is False')
parser.add_argument('--keepH-rmsd', action='store_true', help='Keep H for RMSD calculation. Default is False')


args = parser.parse_args()

if args.keepH_rmsd:
    removeH_rmsd = False
else:
    removeH_rmsd = True

log = open('align.log', mode='w', buffering=1)

ref_elements, ref_coordinates = xyzutils.read_xyz_file(args.reference)
ref_mol = rdkitutils.convert_coordinates_to_mols(ref_elements, ref_coordinates)

ensemble = xyzutils.read_xyz_ensemble(args.ensemble)
ens_elements = [_['elements'] for _ in ensemble.values()]
ens_coordinates = [_['coordinates'] for _ in ensemble.values()]
ens_header = [_['header'] for _ in ensemble.values()]
ens_mols = [rdkitutils.convert_coordinates_to_mols(ele, coords, removeHs=args.removeH_struc) for ele, coords in zip(ens_elements, ens_coordinates)]

atomMap = rdkitutils.get_maximum_substructure_matches(ens_mols, max_matches=1000000, removeHs=removeH_rmsd)
log.write(f'Number of atom maps: {len(atomMap)}\n')

log.write(f'Aligning ensemble {args.ensemble} to reference {args.reference}\n')
calculated_rmsd = {}
for i, probe_mol in enumerate(ens_mols):

    rmsd = rdkitutils.rmsd(probe_mol, ref_mol, atomMap)
    calculated_rmsd[f'{ens_header[i]}'] = rmsd

log.write('Writing aligned ensemble to file\n')
with open(f"{args.ensemble.rsplit('.')[0]}_aligned.xyz", mode='w') as f:
    for mol in ens_mols:
        f.write(Chem.MolToXYZBlock(mol))

calculated_rmsd = pd.Series(calculated_rmsd, name='rmsd', dtype=float)

log.write('Writing calculated RMSD to file\n')
calculated_rmsd.to_csv('rmsd.csv', float_format='%.5f')

log.write('Finished\n\n')
log.close()