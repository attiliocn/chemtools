#!/usr/bin/env python3

import argparse
from modules import geometry, xyzutils, rdkitutils
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('ensemble', help='XYZ Ensemble File')
args = parser.parse_args()

basename, extension = args.ensemble.rsplit('.', 1)

ensemble = xyzutils.read_xyz_ensemble(args.ensemble)
numconfs = len(ensemble)
ens_elements = [_['elements'] for _ in ensemble.values()]
ens_coordinates = [_['coordinates'] for _ in ensemble.values()]
ens_header = [_['header'] for _ in ensemble.values()]
ens_mols = [rdkitutils.convert_coordinates_to_mols(ele, coords) for ele, coords in zip(ens_elements, ens_coordinates)]

rmsd_distance_matrix = rdkitutils.rmsd_matrix_parallel(ens_mols)
np.savetxt(f"{basename}.csv", rmsd_distance_matrix, delimiter=",")