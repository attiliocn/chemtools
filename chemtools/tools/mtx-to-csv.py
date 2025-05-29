#!/usr/bin/env python3

import argparse
import numpy as np
from scipy.io import mmread

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--files', nargs='+', help='.mtx files')
    args = parser.parse_args()

    for file in args.files:
        sparse_matrix = mmread(file)
        dense_matrix = sparse_matrix.toarray()

        output_filename = file.replace('.mtx', '.csv')
        np.savetxt(output_filename, dense_matrix, delimiter=',', fmt='%.3f',)
