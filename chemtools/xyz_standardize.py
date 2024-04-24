#!/usr/bin/env python3

import numpy as np
import sys

from util.xyzutils import read_xyz_file, write_xyz_file

for file in sys.argv[1:]:
    elements, coordinates = read_xyz_file(file)
    coordinates = coordinates.round(4)
    write_xyz_file(elements, coordinates, file)

