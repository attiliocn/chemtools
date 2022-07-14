#!/usr/bin/env python3

from argparse import ArgumentParser
import numpy as np

argparser = ArgumentParser()
argparser.add_argument('files', help='xyz file to mirror', nargs='+')

args = argparser.parse_args()

def read_xyz(fpath):
    xyz_file = open(fpath, mode='r')

    xyz_file.readline() # skip header1 (number of atoms)
    xyz_file.readline() # skip header2 (filename)
    xyz_content = [i.strip().split() for i in xyz_file.readlines()]
    xyz_file.close()

    xyz_array = np.array(xyz_content).reshape(-1,4)
    
    elements = xyz_array[:,0].reshape(-1,1)
    coords = xyz_array[:,1:].astype(np.float64)

    return elements,coords

def mirror_xyz(coords):
    return coords * -1

def write_xyz(elements, coords, filename):
    with open(filename, mode='w') as f:
        f.write(f"{len(coords)}\n")
        f.write(f"{filename}\n")
        for i,j in zip(elements, coords):
            coordinate_line = f'{str(i[0])} {" ".join([str(i) for i in j])}'
            f.write("{: >3} {: >10} {: >10} {: >10}\n".format(*coordinate_line.split()))

    
for file in args.files:
    basename = file.split("/")[-1].split(".")[0]
    elements, coords = read_xyz(file)
    coords_mirror = mirror_xyz(coords)
    write_xyz(elements, coords_mirror, f"{basename}_mirror.xyz")