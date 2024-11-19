#!/usr/bin/env python3
import numpy as np
import argparse
import pandas as pd

from modules.geometry import measure_distance, measure_angle, measure_dihedral_angle
from modules.xyzutils import read_xyz_file

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--files', nargs='+')
parser.add_argument('-m', '--measurement', nargs='+',type=int,action='append')
parser.add_argument('--ignore-atom', action='store_true')
args = parser.parse_args()

def create_measurement_tag(atoms, elements, ignore_atom=False):
    tag = []
    for atom in atoms:
        if not ignore_atom:
            tag.append(f"{elements[atom-1][0]}{atom}")
        elif ignore_atom:
            tag.append(f"{atom} ")
    return "".join(tag)

measurements = {}
for file in args.files:
    measurements[file] = {}
    elements, coordinates = read_xyz_file(file)
    for measurement in args.measurement:
        tag = create_measurement_tag(measurement,elements,args.ignore_atom)
        point_coordinates = [coordinates[i-1] for i in measurement]
        if len(measurement) == 2:
            measurements[file][f"B {tag}"] = measure_distance(*point_coordinates)
        elif len(measurement) == 3:
            measurements[file][f"A {tag}"] = measure_angle(*point_coordinates)
        elif len(measurement) == 4:
            measurements[file][f"D {tag}"] = measure_dihedral_angle(*point_coordinates)

df = pd.DataFrame(measurements).T
df.to_csv('measurements.tsv', sep='\t')