#!/usr/bin/env python3

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-f','--files',nargs='+',required=True)
parser.add_argument('-s','--swap',nargs='+',type=int,action='append',required=True)
args = parser.parse_args()

all_numbers = []
for element in args.swap:
    if isinstance(element,list):
        for i in element:
            all_numbers.append(i-1)

src_numbering = all_numbers[0::2]
dest_numbering = all_numbers[1::2]


def reorder_xyz(xyz_matrix, src_numbering, dest_numbering):
    n_atoms = len(xyz_matrix)
    modified_matrix = np.zeros(n_atoms*4).reshape(n_atoms,4)
    modified_matrix = modified_matrix.astype(xyz_matrix.dtype)

    unmodified_in_source = []
    for i in range(n_atoms):
        if i not in src_numbering:
            unmodified_in_source.append(i)
    unmodified_in_source = np.array(unmodified_in_source)
    print(f"SRC: {np.array(src_numbering)+1}")
    print(unmodified_in_source+1)
    print(len(unmodified_in_source))

    empty_in_destination = []
    for i in range(n_atoms):
        if i not in dest_numbering:
            empty_in_destination.append(i)
    empty_in_destination = np.array(empty_in_destination)
    
    print(f"DEST: {np.array(dest_numbering)+1}")
    print(empty_in_destination+1)
    print(len(empty_in_destination))

    for src, dest in zip(empty_in_destination, unmodified_in_source):
        print(f"S:{src+1} --> D:{dest+1}")
        modified_matrix[src] = xyz_matrix[dest]

    for src, dest in zip(src_numbering, dest_numbering):
       print(f"S:{src+1} --> D:{dest+1}")
       modified_matrix[dest] = xyz_matrix[src]
    
    return modified_matrix

def reorder_matrix_rows(files):
    for file in files:
        with open(file) as f:
            content = f.readlines()
            
        xyz_header = [i.strip() for i in content[0:2]]
        xyz_matrix = [i.split() for i in content[2:]]
        xyz_matrix = np.array(xyz_matrix)
        xyz_matrix = reorder_xyz(xyz_matrix, src_numbering, dest_numbering)

        with open(f"{file.replace('.xyz', '_numbered.xyz')}", mode='w') as f:
            f.write("\n".join(xyz_header))
            f.write("\n")
            for i in xyz_matrix:
                f.write("{: >3} {: >20} {: >20} {: >20}\n".format(*i))

if __name__ == '__main__':
    reorder_matrix_rows(args.files)