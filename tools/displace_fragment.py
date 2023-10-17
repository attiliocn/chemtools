#!/usr/bin/env python3

import numpy as np
import argparse
from morfeus.geometry import kabsch_rotation_matrix

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', help='XMOL .xyz file')
parser.add_argument('--fragment', help='Fragment atom numbers (1-indexed)', type=str, required=True)
parser.add_argument('--stepsize', help='Displacement size (A)', type=float, required=True)
parser.add_argument('--steps', help='Number of displacement steps', type=int, required=True)
parser.add_argument('--dual', help="Use +step and -step sizes", action='store_true')
parser.add_argument('--axis', help='Use two atoms to define the translation axis (Z). 1-indexed', type=int, nargs=2,)

args = parser.parse_args()

def expand_slice_from_string(str_slices):
    parsed_slice = []
    intervals = str_slices.split(',')

    for step in intervals:
        try:
            parsed_slice.append(int(step))
        except:
            range_start, range_stop = [int(i) for i in step.split('-')]
            parsed_slice.extend(range(range_start,range_stop+1))

    parsed_slice = np.array(parsed_slice)
    parsed_slice.sort()
    return parsed_slice

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

def build_xyz_content(elements, coordinates):
    xyz_content = []
    xyz_content.append(f"{str(len(coordinates))}\n")
    xyz_content.append('\n')
    for i in np.concatenate((elements, coordinates), axis=1):
        xyz_content.append("{: >3} {: >15} {: >15} {: >15}\n".format(*i))
    xyz_content = ''.join(xyz_content)
    return xyz_content


elements, coordinates = read_xyz_file(args.file)

if args.axis:
    z_axis = np.array([[0.0, 0.0, 1.0]])
    origin = coordinates[args.axis[0]-1]
    coordinates = coordinates - origin
    axis_bond = coordinates[args.axis[1]-1] - coordinates[args.axis[0]-1]
    R = kabsch_rotation_matrix(axis_bond.reshape(1,-1), z_axis.reshape(1,-1), center=False)
    coordinates = (R @ coordinates.T).T
    displacement = np.array([0.,0.,args.stepsize])
else:
    displacement = args.stepsize

fragment_atoms = expand_slice_from_string(args.fragment) - 1
translate_steps = np.array(range(args.steps+1))[1:]
translate_stepsizes = np.array([i * displacement for i in translate_steps])

if args.dual:
    translate_stepsizes = np.concatenate(
        [
            (translate_stepsizes[::-1]*-1).flatten(),
            [0,0,0],
            (translate_stepsizes).flatten()
        ]
    )
    translate_stepsizes = translate_stepsizes.reshape(-1,3)
else:
    translate_stepsizes = np.concatenate(
        [
            [0,0,0],
            (translate_stepsizes).flatten()
        ]
    )
    translate_stepsizes = translate_stepsizes.reshape(-1,3)

translate_map = [(i+1,j.round(3)) for i,j in enumerate(translate_stepsizes)]

for step, stepsize in translate_map:
    current_coordinates = coordinates.copy()
    current_coordinates[fragment_atoms,:] = current_coordinates[fragment_atoms,:] + stepsize
    current_xmol = build_xyz_content(elements, current_coordinates.round(5))
    with open(args.file.replace('.xyz', f"_scan_{step}.xyz"), mode='w') as f:
        f.write(current_xmol)