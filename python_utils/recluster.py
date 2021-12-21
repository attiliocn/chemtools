#!/usr/bin/env python

import os
from natsort import natsorted

xyz_files = list()
xyz_files_bases = list()
with os.scandir('.') as dir_content:
    for content in dir_content:
        if content.name.endswith('.xyz'):
            xyz_files.append(content.name)
            xyz_base = content.name.rsplit('_', 1)[0]
            xyz_files_bases.append(xyz_base)


xyz_files_bases = natsorted(list(set(xyz_files_bases)))
for conf_family in xyz_files_bases:
    conformers = [i for i in xyz_files if conf_family in i]
    conformers = natsorted(conformers)

    with open(f"{conf_family}_confs.xyz", 'w') as f:
        for conformer in conformers:
            conformer_content = open(conformer)
            f.write(conformer_content.read())
            conformer_content.close()
            

    
