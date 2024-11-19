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

# update destination with previous modified rows
for i,j in enumerate(src_numbering):
    for k,l in enumerate(dest_numbering):
        if k > l and l == j:
            dest_numbering[k] = dest_numbering[i]

def reorder_matrix_rows(files):
    for file in files:
        with open(file) as f:
            content = f.readlines()
            xyz_header = [i.strip() for i in content[0:2]]
            xyz_matrix = [i.split() for i in content[2:]]

        xyz_matrix = np.array(xyz_matrix)
        for swap_pair in zip(src_numbering, dest_numbering):
            xyz_matrix[list(swap_pair)] = xyz_matrix[list(swap_pair)[::-1]]

        with open(f"{file.replace('.xyz', '_numbered.xyz')}", mode='w') as f:
            f.write("\n".join(xyz_header))
            f.write("\n")
            for i in xyz_matrix:
                f.write("{: >3} {: >20} {: >20} {: >20}\n".format(*i))

if __name__ == '__main__':
    reorder_matrix_rows(args.files)