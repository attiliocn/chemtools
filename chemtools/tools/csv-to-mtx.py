#!/usr/bin/env python3

import argparse

import numpy as np
from scipy.io import mmread, mmwrite
from scipy.sparse import coo_matrix

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--files', nargs='+', help='.csv files')
    args = parser.parse_args()

    for file in args.files:
        dense_matrix = np.loadtxt(file, delimiter=',')
        compressed_mtx = coo_matrix(np.tril(dense_matrix))

        output_filename = file.replace('.csv', '.mtx')
        mmwrite(f"{output_filename}", compressed_mtx)