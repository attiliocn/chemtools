#!/usr/bin/env python3

import numpy as np
from scipy.sparse import coo_matrix
from scipy.io import mmwrite
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('files', nargs='+', help='Symmetric lower triangular matrix in csv format')
args = parser.parse_args()

for file in args.files:
    basename, extension = file.rsplit('.', 1)
    original_mtx = np.genfromtxt(file, delimiter=',')
    compressed_mtx = coo_matrix(np.tril(original_mtx.round(3)))
    mmwrite(f"{basename}.mtx", compressed_mtx)

