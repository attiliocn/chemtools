#!/usr/bin/env python
from rdkit import Chem
from openbabel import openbabel
import sys
import argparse
import numpy as np

argparser = argparse.ArgumentParser()
argparser.add_argument('-f', '--file', type=str)
argparser.add_argument('-s', '--smarts', type=str)
argparser.add_argument('-n', '--numbers', nargs='+', type=str)
args = argparser.parse_args()

def obabel(input_string, input_format, output_format):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(input_format, output_format)
    obabel_mol = openbabel.OBMol()
    obConversion.ReadString(obabel_mol,input_string)
    return obConversion.WriteString(obabel_mol)

def get_matches_fragment(xyz_path, fragment):
    xyz_file = open(xyz_path)
    xyz_content = xyz_file.readlines()
    xyz_content = "".join(xyz_content)
    xyz_file.close()
        
    mol2_content = obabel(xyz_content, 'xyz', 'mol2')
    
    mol = Chem.rdmolfiles.MolFromMol2Block(mol2_content, removeHs=False, sanitize=False)
    smarts = Chem.MolFromSmarts(fragment)
    matches = mol.GetSubstructMatches(smarts)   
    return matches

smarts = args.smarts
numbers = args.numbers
file = args.file

matches = get_matches_fragment(file, smarts)
matches = np.array(matches)+1

if matches.any():
    matches = matches[0]
else:
    print(f"Error: SMARTS {args.smarts} does not match for file {args.file}")


if len(matches) == len(numbers):
    exchange_str = ''
    for exchange in zip(matches, numbers):
        if exchange[1] == '_':
            continue
        else:
            exchange_str = exchange_str + "-s {} {} ".format(exchange[0], exchange[1])

    print(f"matrix_n2.py -f {args.file} {exchange_str}")
else:
    print(f'Error: Length mismatch: len(MATCHES) = {len(matches)} and len(NUMBERS) = {len(numbers)}')