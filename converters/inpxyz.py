#!/usr/bin/env python3

import numpy as np
import re
import sys

for file in sys.argv[1:]:
    with open(file) as f:
        content = f.readlines()

    content = [i.strip() for i in content]
    xmol_data = []
    for line in content:
        xyz_pattern = re.search(r"[\s\t]?([a-zA-Z]{1,2})[\s\t]+(-?[0-9]+\.[0-9]+)[\s\t]+(-?[0-9]+\.[0-9]+)[\s\t]+(-?[0-9]+\.[0-9]+)", line) 
        if xyz_pattern:
            xmol_data.append(line)

    with open(f"{file.rsplit('.',1)[0]}.xyz", mode='w') as f: 
        f.write(f"{len(xmol_data)}\n")
        f.write(f"{file}\n")
        for i in xmol_data:
            f.write(f"{i}\n")
    

