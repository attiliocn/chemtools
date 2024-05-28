#!/usr/bin/env python3

import argparse
from modules import geometry, xyzutils
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', help='XYZ Ensemble File')
args = parser.parse_args()

ensemble = xyzutils.read_xyz_ensemble(args.file)
coordinates_all = [_['coordinates'] for _ in ensemble.values()]

rmsd_distance_matrix = geometry.rmsd_matrix_parallel(coordinates_all)
np.savetxt("rmsd.csv", rmsd_distance_matrix, delimiter=",")