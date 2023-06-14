#!/usr/bin/env python3

import numpy as np 
import os
from chemparser.parser_gaussian16 import GaussianOutput

def build_xyz_content(elements, coordinates, title=""):
    xyz_content = []
    xyz_content.append(f"{str(len(coordinates))}\n")
    xyz_content.append(f'{title}\n')
    for i in np.concatenate((elements, coordinates), axis=1):
        xyz_content.append("{: >3} {: >10} {: >10} {: >10}\n".format(*i))
    xyz_content = ''.join(xyz_content)
    return xyz_content

with os.scandir() as d:
    for f in d:
        if f.name.endswith('log'):
            gaussian_output = GaussianOutput(f.name)
            geometries = gaussian_output.get_geometries()

            elements, coordinates = geometries[0][:,0], geometries[0][:,1:]
            elements = elements.reshape(-1,1)
            geom0 = build_xyz_content(elements, coordinates)

            elements, coordinates = geometries[-1][:,0], geometries[-1][:,1:]
            elements = elements.reshape(-1,1)
            geom1 = build_xyz_content(elements, coordinates)

            with open(f"{f.name.replace('.log', '_0.xyz')}", mode='w') as output:
                output.write(geom0)

            with open(f"{f.name.replace('.log', '_1.xyz')}", mode='w') as output:
                output.write(geom1)
