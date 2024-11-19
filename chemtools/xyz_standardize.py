#!/usr/bin/env python3

import numpy as np
import sys

from modules.xyzutils import read_xyz_file, build_xyz_file

for file in sys.argv[1:]:
    elements, coordinates = read_xyz_file(file)
    coordinates = coordinates.round(10)

    _ = build_xyz_file(elements, coordinates, header=file)
    with open(file, mode='w') as f:
        f.write(_)