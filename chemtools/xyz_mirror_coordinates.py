#!/usr/bin/env python3

from argparse import ArgumentParser
import numpy as np
from util.xyzutils import read_xyz_file, write_xyz_file

argparser = ArgumentParser()
argparser.add_argument('files', help='xyz file to mirror', nargs='+')
args = argparser.parse_args()

def mirror_xyz(coords):
    return coords * -1
  
for file in args.files:
    basename = file.split("/")[-1].split(".")[0]
    elements, coords = read_xyz_file(file)
    coords_mirror = mirror_xyz(coords)
    write_xyz_file(elements, coords_mirror, f"{basename}_mirror.xyz")