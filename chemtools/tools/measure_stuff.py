#!/usr/bin/env python3
import numpy as np
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--files', nargs='+')
parser.add_argument('-m', '--measurement', nargs='+',type=int,action='append')
parser.add_argument('--ignore-atom', action='store_true')

args = parser.parse_args()

def measure_distance(p0,p1):
    b0 = p1 - p0
    distance = np.linalg.norm(b0)
    return distance

def measure_angle(p0,p1,p2):
    v1 = p0 - p1
    v2 = p2 - p1
    inner = np.inner(v1, v2)
    norms = np.linalg.norm(v1) * np.linalg.norm(v2)
    cos = inner / norms
    rad = np.arccos(np.clip(cos, -1.0, 1.0))
    deg = np.rad2deg(rad)
    return deg

def measure_dihedral_angle(p0,p1,p2,p3):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

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