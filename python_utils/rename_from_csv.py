#!/usr/bin/env python3

import os

with open('new_names.csv') as f:
    for line in f:
        data = line.strip()
        old_name = data.split(',',1)[0]
        new_name = data.split(',',1)[1]
        if os.path.isfile(old_name):
            print(f"Rename {old_name} to {new_name}")
            os.rename(old_name,new_name)
