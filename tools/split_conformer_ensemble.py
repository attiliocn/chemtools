#!/usr/bin/env python3

import argparse, os, re

argparser = argparse.ArgumentParser()
argparser.add_argument('files', nargs='+')
argparser.add_argument('--keep-name', action='store_true')
argparser.add_argument('--split-char', default='_')
argparser.add_argument('--split-level', default=1, type=int)
args = argparser.parse_args()

for ensemble_filename in args.files:
    if args.keep_name:
        with open(ensemble_filename) as f:
            ensemble_content = f.readlines()

        num_atoms = int(ensemble_content[0])
        ensemble_split = []
        current_line = 0       
        while current_line < len(ensemble_content):
            ensemble_split.append(ensemble_content[current_line:(current_line+num_atoms+2)])
            current_line += (num_atoms + 2)
        for molecule in ensemble_split:
            filename = f"{molecule[1].split('.')[0]}.xyz"
            with open(filename, mode='w') as f:
                f.write("".join(molecule))
    else:
        basename = ensemble_filename.rsplit(args.split_char, args.split_level)[0]
        output_filename = f"{basename}_.xyz"
        os.system(f'obabel -ixyz {ensemble_filename} -oxyz -m -O {output_filename}')