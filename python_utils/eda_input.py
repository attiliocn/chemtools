#!/usr/bin/env python
import argparse
from shutil import move


parser = argparse.ArgumentParser()
parser.add_argument('xyz_files', nargs='+')
parser.add_argument('-frag1', nargs='+', required=True)
parser.add_argument('--in-place', action='store_true')
args = parser.parse_args()
frag1_coords = [int(i)-1 for i in args.frag1]

for xyz_file in args.xyz_files:
    with open('temp.xyz', mode='w') as new_xyz:
        with open(xyz_file) as f:
            new_xyz.write(f.readline()) # write number of atoms
            new_xyz.write(f.readline()) # write XYZ file title
            input_content = f.readlines()
            for i, coord in enumerate(input_content):
                coord = coord.split()
                if i in frag1_coords:
                    coord[0] = f"{coord[0]}(fragment=1)"
                    coord_string = '    '.join(coord)
                    new_xyz.write(f"{coord_string}\n")
                else:
                    coord[0] = f"{coord[0]}(fragment=2)"
                    coord_string = '    '.join(coord)
                    new_xyz.write(f"{coord_string}\n")
    if args.in_place:
        move('temp.xyz', xyz_file)
    else:
        move('temp.xyz', f"mod_{xyz_file}")