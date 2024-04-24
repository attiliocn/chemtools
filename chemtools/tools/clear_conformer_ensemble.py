#!/usr/bin/env python3

import os
import argparse

import numpy as np
import pandas as pd
from itertools import combinations

from spyrmsd import rmsd
from parse_xyz import parse_xyz_ensemble as read_xyz
from parse_xyz import write_xyz_file as write_xyz

argparser = argparse.ArgumentParser()
argparser.add_argument('files', help='Ensemble File in XYZ format', nargs='+')
argparser.add_argument('--write', help='Write the unclustered ensemble to a XYZ file', action='store_true')
args = argparser.parse_args()

def calculateDistanceMatrix(conf_ensemble):
    matrixSize = len(conf_ensemble)
    rmsdMatrix = np.zeros(matrixSize**2).reshape(matrixSize,matrixSize)
    
    comb = combinations(range(matrixSize),2)
    comb = list(comb)
    
    for matrixElement in comb:
        i = matrixElement[0]
        j = matrixElement[1]
      
        rmsd_entry = rmsd.hrmsd(
            conf_ensemble[i]['coordinates'], conf_ensemble[j]['coordinates'],
            conf_ensemble[i]['atomic_numbers'], conf_ensemble[j]['atomic_numbers'],
        )
        
        rmsdMatrix[i][j] = rmsd_entry
        rmsdMatrix[j][i] = rmsd_entry
        
    return rmsdMatrix

liner = '#'*80
dropped_indexes = dict()

if os.path.exists('drops.log'): os.remove('drops.log')

for ensemble in args.files:
    ens_base = ensemble.split('.')[0]
    
    #load molecules
    molecules = read_xyz(ensemble)
    matrix = calculateDistanceMatrix(molecules)
    matrix = pd.DataFrame(matrix)
    matrix.to_csv('drops.csv')
    
    #process the unclustered distance matrix based on an RMSD threshold
    with open('drops.log', mode='a') as f:
        f.write(f'{liner}\n')
        f.write(f'Processing {ens_base}\n')
        f.write(f'Number of conformers {len(molecules)}\n')
        
        dump = matrix.copy()
        droppeds = list() #0 indexed list of conformers to drop
        for conf_id in matrix.index.values:

            if conf_id in droppeds:
                f.write(f"ID {conf_id+1} already dropped\n")
            else:
                to_drop = matrix[matrix[conf_id] < 0.30].index
                to_drop = list(to_drop)

                to_drop.remove(conf_id)

                for dropped in droppeds:
                    if dropped in to_drop:
                        to_drop.remove(dropped)
                        
                f.write(f"Current conf: {conf_id+1}, will drop the confs: {[x+1 for x in to_drop]}\n")
                for conf_to_drop in to_drop:
                    droppeds.append(conf_to_drop)
                    dump.drop(conf_to_drop, axis=0, inplace=True)
                    dump.drop(conf_to_drop, axis=1, inplace=True)
        
        dropped_indexes[ens_base] = [x+1 for x in droppeds] #1 indexed list of dropped conformers
        
        f.write(f"{str([x+1 for x in droppeds])}\n")
        f.write(f'Number of conformers {len(molecules)}\n')
        f.write(f"Dropped {len(droppeds)}\n")
        f.write(f'Final conformers {len(molecules) - len(droppeds)}\n')
        f.write(f"Processing of {ens_base} is done\n")
        f.write(f'{liner}\n\n')

    mols_ids = set(molecules.keys())
    mols_drop_ids = set(droppeds)
    remaining_mols = list(mols_ids.difference(mols_drop_ids))

    output_ensemble = f"{ens_base}_drop.xyz"
    output_ensemble_exists = os.path.isfile(output_ensemble)

    if output_ensemble_exists: os.remove(output_ensemble)

    for id in remaining_mols:
        write_xyz(
            molecules[id]['elements'],
            molecules[id]['coordinates'],
            molecules[id]['header'],
            filename=output_ensemble
        )

        
