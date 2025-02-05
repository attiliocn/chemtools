#!/usr/bin/env python3

import argparse 
from modules import xyzutils
import os

parser = argparse.ArgumentParser()
parser.add_argument('input', help='3-columns .csv input file. TS Geometry, TS Energy, FW Role. Must include header line')
args = parser.parse_args()

with open(args.input) as f:
    input_file_content = [l.strip().split(',',3) for l in f]

ts_recipes = []
for content in input_file_content[1:]:
    ts_recipes.append(
        (
            content[0], 
            float(content[1]),
            int(content[2])
        )
    )

for recipe in ts_recipes:
    basename = recipe[0].replace('.xyz','')

    ts_geometry = xyzutils.read_xyz_ensemble(recipe[0])
    ts_geometry[0]['stringContent'] = xyzutils.build_xyz_file(
        ts_geometry[0]['elements'],
        ts_geometry[0]['coordinates'],
        header=f'TS Geometry E {recipe[1]}'
    )

    rv_trajectory = xyzutils.read_xyz_ensemble(f'{basename}_rv.trj')

    if os.path.exists(f'downhill/{basename}_rv_down.trj'):
        rv_downhill_trajectory = xyzutils.read_xyz_ensemble(f'downhill/{basename}_rv_down.trj')
    else:
        rv_downhill_trajectory = None

    if os.path.exists(f'endpoints/{basename}_rv.xyz'):
        rv_endpoint_geometry = xyzutils.read_xyz_ensemble(f'endpoints/{basename}_rv.xyz')
        # rv_endpoint_geometry[0]['stringContent'] = xyzutils.build_xyz_file(
        #     rv_endpoint_geometry[0]['elements'],
        #     rv_endpoint_geometry[0]['coordinates'],
        #     header=f'RV Endpoint REPLACE_ENERGY'
        # )
    else:
        rv_endpoint_geometry = None

    fw_trajectory = xyzutils.read_xyz_ensemble(f'{basename}_fw.trj')

    if os.path.exists(f'downhill/{basename}_fw_down.trj'):
        fw_downhill_trajectory = xyzutils.read_xyz_ensemble(f'downhill/{basename}_fw_down.trj')
    else:
        fw_downhill_trajectory = None

    if os.path.exists(f'endpoints/{basename}_rv.xyz'):
        fw_endpoint_geometry = xyzutils.read_xyz_ensemble(f'endpoints/{basename}_fw.xyz')
        # fw_endpoint_geometry[0]['stringContent'] = xyzutils.build_xyz_file(
        #     fw_endpoint_geometry[0]['elements'],
        #     fw_endpoint_geometry[0]['coordinates'],
        #     header=f'FW Endpoint REPLACE_ENERGY'
        # )
    else:
        fw_endpoint_geometry = None

    if recipe[2] == 0:
        # then fw is the trajectory toward reactants
        # fw trajectory must be reversed
        # the full path is fw(reversed) -> ts -> rv(unmodified)
        with open(f'{basename}.irc', mode='w') as f:
            if fw_endpoint_geometry:
                f.write(fw_endpoint_geometry[0]['stringContent'])

            if fw_downhill_trajectory:
                for mol_id in list(fw_downhill_trajectory.keys())[:0:-1]:
                    f.write(fw_downhill_trajectory[mol_id]['stringContent'])

            for mol_id in list(fw_trajectory.keys())[::-1]:
                f.write(fw_trajectory[mol_id]['stringContent'])
            
            f.write(ts_geometry[0]['stringContent'])
            
            for mol_id in list(rv_trajectory.keys()):
                f.write(rv_trajectory[mol_id]['stringContent'])

            if rv_downhill_trajectory:
                for mol_id in list(rv_downhill_trajectory.keys())[1:]:
                    f.write(rv_downhill_trajectory[mol_id]['stringContent'])
            
            if rv_endpoint_geometry:
                f.write(rv_endpoint_geometry[0]['stringContent'])

    elif recipe[2] == 1:
        # then fw is the trajectory toward products
        # fw trajectory remains unmodified
        # rv trajectory must be reversed
        # the full path is then rv(reversed) -> ts -> fw(unmodified)
        with open(f'{basename}.irc', mode='w') as f:
            if rv_endpoint_geometry:
                f.write(rv_endpoint_geometry[0]['stringContent'])

            if rv_downhill_trajectory:
                for mol_id in list(rv_downhill_trajectory.keys())[:0:-1]:
                    f.write(rv_downhill_trajectory[mol_id]['stringContent'])

            for mol_id in list(rv_trajectory.keys())[::-1]:
                f.write(rv_trajectory[mol_id]['stringContent'])
            
            f.write(ts_geometry[0]['stringContent'])

            for mol_id in list(fw_trajectory.keys()):
                f.write(fw_trajectory[mol_id]['stringContent'])

            if fw_downhill_trajectory:
                for mol_id in list(fw_downhill_trajectory.keys())[1:]:
                    f.write(fw_downhill_trajectory[mol_id]['stringContent'])

            if fw_endpoint_geometry:
                f.write(fw_endpoint_geometry[0]['stringContent'])
            
