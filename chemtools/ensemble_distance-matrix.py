#!/usr/bin/env python3

import argparse
from modules import geometry, xyzutils
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('ensemble', help='XYZ Ensemble File')
parser.add_argument('--all-atoms', action='store_true', help='Include H atoms in the RMSD calculation. Default is use only heavy atoms')
args = parser.parse_args()

basename, extension = args.ensemble.rsplit('.', 1)

ensemble = xyzutils.read_xyz_ensemble(args.ensemble)
coordinates_all = [_['coordinates'] for _ in ensemble.values()]
elements_all = [_['elements'] for _ in ensemble.values()]

if args.all_atoms:
        coordinates_rmsd = coordinates_all
else:
    heavy_atoms = np.where(elements_all[0] != 'H')[0]
    coordinates_rmsd = [i[heavy_atoms,:] for i in coordinates_all]

rmsd_distance_matrix = geometry.rmsd_matrix_parallel(coordinates_rmsd)
np.savetxt(f"{basename}.csv", rmsd_distance_matrix, delimiter=",")