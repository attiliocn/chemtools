#!/usr/bin/env python3

import argparse

from modules.xyzutils import read_xyz_ensemble, build_xyz_file

parser = argparse.ArgumentParser()
parser.add_argument('file', nargs='+', help='XYZ Files')
parser.add_argument('--geometry', type=int, default=1, help='0-indexed geometry number. Also, -1 extract the last geometry and so on')
parser.add_argument('--prefix', type=int, help='Remove n prefixes from filename (separated by underscore)')
args = parser.parse_args()

def main(filename):

    if args.prefix:
        new_filename = filename.split('_')
        keep_indexes = len(new_filename) - args.prefix
        new_filename = "_".join(new_filename[:keep_indexes])
    else:
        new_filename = filename.split('.')[0]
    new_filename = f"{new_filename}_{args.geometry}.xyz"

    ensemble = read_xyz_ensemble(filename)

    molecules_idx = list(ensemble.keys())
    request = molecules_idx[args.geometry]

    molecule = ensemble[request]
    molecule['coordinates'] = molecule['coordinates'].round(5)
    _ = build_xyz_file(
        molecule['elements'], 
        molecule['coordinates'],
        header=f"{filename} geometry no {args.geometry}"
    )

    with open(new_filename, mode='w') as f:
        f.write(_)

for file in args.file:
    main(file)