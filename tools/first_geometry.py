#!/usr/bin/env python3

import argparse

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

    with open(filename) as f:
        all_geometries = f.readlines()

    num_atoms = int(all_geometries[0].strip())

    start_line = (num_atoms + 2) * (args.geometry - 1) 
    end_line = (num_atoms + 2) * (args.geometry)

    print("Extracting data from the following lines:")
    print(f"Start: {start_line+1}\nEnd: {end_line+1}")
    print(f"New filename: {new_filename}")

    with open(new_filename, mode='w') as f:
        for i in all_geometries[start_line:end_line]:
            f.write(i)

for file in args.file:
    main(file)
