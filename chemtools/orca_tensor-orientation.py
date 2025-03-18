#! /usr/bin/env python3

import numpy as np
import re
import argparse

def get_tensors_orientations(output:str):
    with open(output) as f:
        content = f.readlines()

    for i, line in enumerate(content):
        if re.search(r'\sOrientation:', line):
            vectors_idx = i+1
            break

    tensor_orientation_raw = content[vectors_idx:vectors_idx+3]

    tensor_orientation = []
    for line in tensor_orientation_raw:
        tensor_orientation.append(re.findall(r'-?\d+\.\d+', line))

    tensor_orientation = np.array(tensor_orientation, dtype=float).T
    return tensor_orientation

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+', help='ORCA 6 output files')
    args = parser.parse_args()

    all_tensors = {}
    for file in args.files:
        all_tensors[file] = get_tensors_orientations(file)

    print(all_tensors)

if __name__ == "__main__":
    main()