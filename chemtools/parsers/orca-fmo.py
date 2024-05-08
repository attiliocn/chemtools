#!/usr/bin/env python3

import argparse
import re
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-f', nargs='+', help='ORCA output files')
parser.add_argument('--units', choices=['Eh', 'eV'], default='Eh')
args = parser.parse_args()

if args.units == 'Eh':
    energy_col = 2
elif args.units == 'eV':
    energy_col = 3

def get_orbital_energies(file):
    orbital_energies = []
    with open(file) as f:
        content = [i.strip() for i in f.readlines()]
        start_lines = []
        for i,l in enumerate(content):
            if re.search('ORBITAL ENERGIES', l):
                start_lines.append(i)

    for sl in start_lines:
        _=[]
        for l in content[sl:]:
            _.append(l)
            if re.search(r'\*\*\*', l):
                break
        _ = [re.sub(r'\s+', ',', i) for i in _]
        _ = [i.split(",") for i in _[4:-2]]
        orbital_energies.append(_)
    
    return orbital_energies        

for output_file in args.f:
    energies = get_orbital_energies(output_file)[-1]
    energies = np.array(energies, dtype=float)

    lumo_id = np.where(energies[:,1] == 0)[0][0]
    homo_id = lumo_id-1

    lumo_energy = energies[lumo_id, energy_col] 
    homo_energy = energies[homo_id, energy_col]
    gap = (lumo_energy - homo_energy).round(4)

    print(f"file,lumo,homo,gap {args.units}")
    print(f"{output_file},{lumo_energy},{homo_energy},{gap}")
