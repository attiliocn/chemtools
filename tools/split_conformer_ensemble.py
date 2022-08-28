#!/usr/bin/env python3

import argparse, os

argparser = argparse.ArgumentParser()
argparser.add_argument('files', nargs='+')
argparser.add_argument('--split-char', default='_')
argparser.add_argument('--split-level', default=1, type=int)
args = argparser.parse_args()

for ensemble_filename in args.files:
    basename = ensemble_filename.rsplit(args.split_char, args.split_level)[0]
    output_filename = f"{basename}_.xyz"
    os.system(f'obabel -ixyz {ensemble_filename} -oxyz -m -O {output_filename}')
