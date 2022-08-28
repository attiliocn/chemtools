#!/usr/bin/env python3

import os, re
from natsort import natsorted
import numpy as np


all_xyz_files = []
with os.scandir('./') as d:
    for f in d:
        if f.name.endswith('xyz'):
            all_xyz_files.append(f.name)

def calculate_relative_energies(energies, factor=1):
    return ((energies-energies.min())*factor)

regex_energy = re.compile('Energy: +(-[0-9]+.[0-9]+)', re.MULTILINE)
all_xyz_groups = list(set([i.rsplit('_',1)[0] for i in all_xyz_files]))
all_xyz_groups.sort()
for group in all_xyz_groups:
    selected_files = natsorted([i for i in all_xyz_files if group in i])
    all_energies = []
    for file in selected_files:
        f = open(file)
        energies = float(regex_energy.findall(f.read())[0])
        all_energies.append(energies)
        f.close()
    
    all_energies = np.array(all_energies)
    relative_energies = calculate_relative_energies(all_energies, factor=1)   
    
    all_energies_parsed = {}
    for i, file in enumerate(selected_files):
        all_energies_parsed[file] = {
            'energy':all_energies[i]/627.5096,
            'rel_energy':relative_energies[i]
        }
    
    for file in all_energies_parsed.keys():
        os.system(f"sed -i '2s/.*/{file} {all_energies_parsed[file]['energy']}/g' {file}")
    