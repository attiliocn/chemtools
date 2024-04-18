#!/usr/bin/env python3

import numpy as np
import sys

def read_xyz_file(filepath):
    with open(filepath) as f:
        f.readline()
        f.readline()
        data = [i.strip().split() for i in f.readlines()]
        elements = [i[0] for i in data]
        coordinates = [i[1:] for i in data]
    elements = np.array(elements).reshape(-1,1)
    coordinates = np.array(coordinates, dtype=float)
    return elements, coordinates

def write_xyz(elements, coords, filename):
    with open(filename, mode='w') as f:
        f.write(f"{len(coords)}\n")
        f.write(f"{filename}\n")
        for i,j in zip(elements, coords):
            coordinate_line = f'{str(i[0])} {" ".join([str(i) for i in j])}'
            f.write("{: >3} {: >10} {: >10} {: >10}\n".format(*coordinate_line.split()))

for file in sys.argv[1:]:
    elements, coordinates = read_xyz_file(file)
    coordinates = coordinates.round(4)
    write_xyz(elements, coordinates, file)

