#!/usr/bin/env python3

import argparse

from chemtools.modules import rdkitutils, xyzutils, kabsch
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description="Align a bond to an axis")
    parser.add_argument('-f', '--file', help='.xyz (xmol) file', required=True)
    parser.add_argument('-a', '--axis', help='The reference axis to align the bond', choices=['x', 'y', 'z'], required=True)
    parser.add_argument('--translate', help='Translate atom1 to origin', action='store_true')
    parser.add_argument('--atom1', help='0 indexed atom number defining the bond origin', nargs=1, type=int)
    parser.add_argument('--atom2', help='0 indexed atom number defining the bond terminus. If n>1 the average cordinate will be used', nargs='+', type=int)
    return parser.parse_args()

def main():
    args = parse_args()
    
    # cartesian base vectors as column-matrices
    x_axis = np.array([1,0,0]).reshape(3,-1)
    y_axis = np.array([0,1,0]).reshape(3,-1)
    z_axis = np.array([0,0,1]).reshape(3,-1)

    cartesian_basis = {
        'x': x_axis,
        'y': y_axis,
        'z': z_axis
    }

    elements, coordinates = xyzutils.read_xyz_file(args.file)

    atom1_coordinates = coordinates[args.atom1]
    atom2_coordinates = coordinates[args.atom2,:]
    atom2_coordinates_average = np.average(atom2_coordinates, axis=0)
   
    bond_vector = atom2_coordinates_average - atom1_coordinates
    bond_vector = bond_vector.reshape(3,-1)

    R,t = kabsch.kabsch_algorithm(
        bond_vector,
        cartesian_basis[args.axis],
        center=False
    )
    coordinates_after_rotation = (R @ coordinates.T).T
    
    if args.translate:
        coordinates_after_rotation -= coordinates_after_rotation[args.atom1]

    final_xyz_file = xyzutils.build_xyz_file(
        elements, coordinates_after_rotation,
        header=f'{args.file} after rotation'
    )
    print(final_xyz_file)

if __name__ == '__main__':
    main()