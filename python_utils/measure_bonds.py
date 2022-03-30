#!/usr/bin/env python
import sys
import numpy as np
import argparse
from termcolor import colored

sys.stdout.reconfigure(line_buffering=True)

parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+')
parser.add_argument('atom', nargs=2, type=int, help='1 indexed atom numbers involved in the desired bond')
parser.add_argument('--verbose', '-v', action='store_true')
args = parser.parse_args()

for file in args.files:
    atoms = list()
    coordinates = list()
    with open(file, mode='r') as f:
        f.readline() # skip number of atoms row
        f.readline() # skip title row
        for line in f:
            atom = line.split()[0]
            atoms.append(atom)
            atom_coordinates = [float(i) for i in line.split()[1:]]
            coordinates.append(atom_coordinates)

    coordinates = np.array(coordinates)
    atom_1 = args.atom[0] - 1 # 0-index atom 1
    atom_2 = args.atom[1] - 1 # 0-index atom 2
    bond_distance = np.linalg.norm(coordinates[atom_1] - coordinates[atom_2])

    if args.verbose:
        print(f"{colored(file, 'blue', attrs=['bold'])}  {bond_distance:.8f}  Å")
    else:
        print(bond_distance)
