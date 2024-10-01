#!/usr/bin/env python3

import argparse

import numpy as np
from modules import geometry, xyzutils

parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+', help='xyz ensemble file (compatible with several files)')
args = parser.parse_args()

def main(filename):
    ensemble = xyzutils.read_xyz_ensemble(filename)
    molecules_idx = list(ensemble.keys())

    with open(f'{filename}.rev', mode='w') as f:
        for idx in molecules_idx[::-1]:
            f.write(ensemble[molecules_idx[idx]]['stringContent'])

for file in args.files:
    main(file)