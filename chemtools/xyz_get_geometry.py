#!/usr/bin/env python3

import argparse

from util.xyzutils import read_xyz_ensemble, build_xyz_content

parser = argparse.ArgumentParser()
parser.add_argument('file', nargs='+', help='XYZ Files')
parser.add_argument('--geometry', type=int, default=1, help='1-indexed geometry number')
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
    molecule = ensemble[args.geometry-1]
    molecule['coordinates'] = molecule['coordinates'].round(5)
    _ = build_xyz_content(
        molecule['elements'], 
        molecule['coordinates'],
        header=f"{args.file} geometry no {args.geometry}"
    )

    with open(new_filename, mode='w') as f:
        f.write(_)

for file in args.file:
    main(file)