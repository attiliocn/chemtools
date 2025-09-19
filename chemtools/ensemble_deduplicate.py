#!/usr/bin/env python3

# TODO: use __main__ to loop through inputs, instead of a fucking 
# _for_ loop

import argparse
import time

import numpy as np
from modules import geometry, xyzutils

from scipy.io import mmread

start_time_global = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+', help='.mtx distance matrix (expected: sparse lower triangular matrix)')
parser.add_argument('--threshold', type=float, default=.25, help='RMSD threshold for duplicate detection. Default is 0.25')
args = parser.parse_args()

log = open('deduplicate.log', mode='w', buffering=1)
log.write(f"THRESHOLD: {args.threshold}\n")

for file in args.files:
    start_time = time.time()

    basename, extension = file.rsplit('.', 1)
    basename_updated = f"{basename}_deduplicated.xyz"
    log.write(f"Current file: {basename}\n")

    # rmsd_distance_matrix = np.genfromtxt(file, delimiter=",")
    rmsd_distance_matrix = mmread(file)
    rmsd_distance_matrix = rmsd_distance_matrix.toarray()
    numconfs = rmsd_distance_matrix.shape[0]
    to_delete, conformer_relationships = geometry.get_duplicates_rmsd_matrix(
        rmsd_distance_matrix, 
        threshold=args.threshold,
        return_relationships=True
    )
    log.write(f"Number of conformers: {numconfs}\n")

    ensemble = xyzutils.read_xyz_ensemble(f"{basename}.xyz")
    ens_elements = [_['elements'] for _ in ensemble.values()]
    ens_coordinates = [_['coordinates'] for _ in ensemble.values()]
    ens_header = [_['header'] for _ in ensemble.values()]
    
    _coordinates = [ens_coordinates[i] for i in range(numconfs) if i not in to_delete]
    _elements = [ens_elements[i] for i in range(numconfs) if i not in to_delete]
    _header = [ens_header[i] for i in range(numconfs) if i not in to_delete]

    log.write(f"Conformers after deduplication: {len(_coordinates)}\n")
    log.write(f"{numconfs - len(_coordinates)} conformers were deleted\n")

    for k,v in conformer_relationships.items():
        log.write(f'Conformer {k+1} will be kept. It is similar to {','.join([str(i+1) for i in v])}\n')

    log.write(f"{len(_coordinates)} conformers will be written to {basename_updated}\n")
    end_time = time.time()
    elapsed_time = end_time - start_time
    log.write(f"Wall clock time: {elapsed_time:.2f} seconds\n\n")

    with open(basename_updated, mode='w') as f:
        for header, elements, coordinates in zip(_header, _elements, _coordinates):  
            _ = xyzutils.build_xyz_file(
                elements, 
                coordinates,
                header=header
            )
            f.write(_)

end_time_global = time.time()
elapsed_time_global = end_time_global - start_time_global
log.write(f"Total wall clock time: {elapsed_time_global:.2f} seconds\n")

log.close()        