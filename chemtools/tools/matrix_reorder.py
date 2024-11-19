#!/usr/bin/env python3

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-f','--files',nargs='+',required=True)
parser.add_argument('-s','--swap',nargs='+',type=int,action='append',required=True)
args = parser.parse_args()

src_numbering = []
for element in args.swap:
    if isinstance(element,list):
        for i in element:
            src_numbering.append(i-1)

dest_numbering = src_numbering.copy()
for i in range(0,len(dest_numbering),2):
    dest_numbering[i], dest_numbering[i+1] = dest_numbering[i+1], dest_numbering[i]

def reorder_matrix_rows(files):
    for file in files:
        with open(file) as f:
            content = f.readlines()
            xyz_header = [i.strip() for i in content[0:2]]
            xyz_matrix = [i.split() for i in content[2:]]

        xyz_matrix = np.array(xyz_matrix)
        xyz_matrix[src_numbering] = xyz_matrix[dest_numbering]

        with open(f"{file.replace('.xyz', '_numbered.xyz')}", mode='w') as f:
            f.write("\n".join(xyz_header))
            f.write("\n")
            for i in xyz_matrix:
                f.write("{: >3} {: >20} {: >20} {: >20}\n".format(*i))

if __name__ == '__main__':
    reorder_matrix_rows(args.files)