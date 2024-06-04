#!/usr/bin/env python3

import argparse
import os
import natsort
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument(
    '-s', 
    '--separator', 
    default='_', 
    help='Delimiting character for sufixes. Default is _'
)
parser.add_argument(
    '-n', 
    '--num-separators', 
    type=int, 
    default=1, 
    help='Number of sufixes to ignore. Default is 1'
)
parser.add_argument(
    '-e', 
    '--file-extension', 
    type=str, 
    help="Count only files with the specified extension. Otherwise ignore the extension"
)
parser.add_argument(
    '--template', 
    help='Text file including the default sufixes'
)
args = parser.parse_args()


files = []
with os.scandir(os.getcwd()) as d:
    for f in d:
        if args.file_extension:
            if not f.name.endswith(args.file_extension):
                continue
        files.append(f.name)
files_parsed = [f.rsplit(args.separator,args.num_separators)[0] for f in files]
files_parsed = natsort.natsorted(files_parsed)
files_counts = Counter(files_parsed)

if args.template:
    with open(args.template) as f:
        to_count = f.readlines()
    to_count = [i.strip() for i in to_count]
    _ = {}
    for entry in to_count:
        try:
            _[entry] = files_counts[entry]
        except KeyError:
            _[entry] = 0   
    files_counts = _

for f, i in files_counts.items():
    print(f"{f:<30}{i:>5}")

