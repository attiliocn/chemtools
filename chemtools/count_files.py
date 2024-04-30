#!/usr/bin/env python3

import argparse
import os
import natsort
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--separator', default='_', help='Delimiting character for sufixes. Default is _')
parser.add_argument('-n', '--num-separators', type=int, default=1, help='Number of sufixes to ignore. Default is 1')
args = parser.parse_args()

files = []
with os.scandir(os.getcwd()) as d:
    for f in d:
        files.append(f.name)

files_parsed = [f.rsplit(args.separator,args.num_separators)[0] for f in files]
files_parsed = natsort.natsorted(files_parsed)

files_counts = Counter(files_parsed)
for f, i in files_counts.items():
    print(f"{f} > {i}")