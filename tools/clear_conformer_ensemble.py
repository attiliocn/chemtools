#!/usr/bin/env python3

import os
import argparse

import numpy as np
import pandas as pd
from itertools import combinations

import matplotlib.pyplot as plt
import seaborn as sns

from spyrmsd import rmsd
from parse_xyz import parse_xyz_ensemble as read_xyz
from parse_xyz import write_xyz_file as write_xyz

argparser = argparse.ArgumentParser()
argparser.add_argument('files', help='Ensemble File in XYZ format', nargs='+')
argparser.add_argument('--fig', help='Create Hierarchical Clustering plots', action='store_true')
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

def save_fig(path, ext='png', dpi=300):
    plt.savefig(f"{path}.{ext}", dpi=dpi)

liner = '#'*80
dropped_indexes = dict()

if os.path.exists('drops.log'): os.remove('drops.log')

for ensemble in args.files:
    ens_base = ensemble.split('.')[0]
    
    #load molecules
    molecules = read_xyz(ensemble)
    matrix = calculateDistanceMatrix(molecules)
    matrix = pd.DataFrame(matrix)
    
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

    if args.fig:
        if not os.path.isdir('fig'): os.mkdir('fig')
        #plot distance matrix in the original order
        f, ax = plt.subplots(figsize=(11, 9))
        sns.heatmap(matrix, cmap='jet')
        save_fig(f'fig/{ens_base}_1-original_order')
        plt.close()
        
        #plot distance matrix after hierarchical clustering
        matrix_clustered = sns.clustermap(matrix, cmap='jet')
        matrix_clustered.savefig(f'fig/{ens_base}_2-cluster_order')
        plt.close()
    
        #plot the distance matrix with hierarchical chlustering after dropping the duplicates 
        matrix_cut = sns.clustermap(dump, cmap='jet')
        matrix_cut.savefig(f'fig/{ens_base}_3-duplicates_removed.png')
        plt.close()

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

        
