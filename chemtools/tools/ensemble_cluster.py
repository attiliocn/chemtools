#!/usr/bin/env python3
import numpy as np
from itertools import combinations
from spyrmsd import rmsd, io

import matplotlib.pyplot as plt
import pandas as pd

import seaborn as sns

def calc_rmsd(n1,n2):
    mol1 = molecules[n1-1]
    mol2 = molecules[n2-1]
    rmsd_calc = rmsd.hrmsd(
            mol1.coordinates, mol2.coordinates,
            mol1.atomicnums, mol2.atomicnums,
        )
    print(rmsd_calc)

def calculateDistanceMatrix(conf_ensemble):
    matrixSize = len(conf_ensemble)
    rmsdMatrix = np.zeros(matrixSize**2).reshape(matrixSize,matrixSize)
    
    comb = combinations(range(matrixSize),2)
    comb = list(comb)
    
    for matrixElement in comb:
        i = matrixElement[0]
        j = matrixElement[1]
        
        rmsd_entry = rmsd.hrmsd(
            conf_ensemble[i].coordinates, conf_ensemble[j].coordinates,
            conf_ensemble[i].atomicnums, conf_ensemble[j].atomicnums,
        )
        
        rmsdMatrix[i][j] = rmsd_entry
        rmsdMatrix[j][i] = rmsd_entry
        
    return rmsdMatrix

def save_fig(path, ext='png', dpi=300):
    plt.savefig(f"{path}.{ext}", dpi=dpi)


confs_filename = 'crest_conformers.xyz'
data_name = 'dump'

molecules = io.loadallmols(f'{confs_filename}')
matrix = calculateDistanceMatrix(molecules)
matrix = pd.DataFrame(matrix)

#matrix_clustered = sns.clustermap(matrix, cmap='jet')
#matrix_clustered.savefig(f'fig/confs-{data_name}_heatmap_clustered.png')


dump = matrix.copy()

droppeds = list()
for conf_id in matrix.index.values:
    
    if conf_id in droppeds:
        print(f"ID {conf_id+1} already dropped")
    else:
        to_drop = matrix[matrix[conf_id] < 0.30].index
        to_drop = list(to_drop)

        to_drop.remove(conf_id)

        for dropped in droppeds:
            if dropped in to_drop:
                to_drop.remove(dropped)

        print(f"Current conf: {conf_id+1}, will drop the confs: {[x+1 for x in to_drop]}")

        for conf_to_drop in to_drop:
            droppeds.append(conf_to_drop)
            dump.drop(conf_to_drop, axis=0, inplace=True)
            dump.drop(conf_to_drop, axis=1, inplace=True)

print([x+1 for x in droppeds])
print(f"Dropped {len(droppeds)}")


#matrix_cut = sns.clustermap(dump, cmap='jet')
#matrix_cut.savefig(f'fig/confs-{data_name}_heatmap_cut.png')
print(len(dump))