#!/usr/bin/env python3

# TODO: make this script operate only in .csv lower triangular distance matrix
# thus avoiding the need to refactor several rmsd engines every time a bug
# comes out

# TODO: use __main__ to loop through inputs, instead of a fucking 
# _for_ loop

import argparse
import time
import os

import numpy as np

from modules import geometry, rdkitutils, xyzutils

start_time_global = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+', help='XYZ Ensemble Files')
parser.add_argument('--threshold', type=float, default=.25, help='RMSD threshold for duplicate detection. Default is 0.25')
args = parser.parse_args()

log = open('deduplicate.log', mode='w', buffering=1)
log.write(f"THRESHOLD: {args.threshold}\n")

for file in args.files:
    basename, extension = file.rsplit('.', 1)
    basename_updated = f"{basename}_deduplicated.{extension}"

    ensemble = xyzutils.read_xyz_ensemble(file)
    numconfs = len(ensemble)
    ens_elements = [_['elements'] for _ in ensemble.values()]
    ens_coordinates = [_['coordinates'] for _ in ensemble.values()]
    ens_header = [_['header'] for _ in ensemble.values()]
    ens_mols = [rdkitutils.convert_coordinates_to_mols(ele, coords) for ele, coords in zip(ens_elements, ens_coordinates)]

    log.write(f"Current file: {basename}\n")
    log.write(f"Number of conformers: {numconfs}\n")

    start_time = time.time()

    log.write("Using parallel RMSD calculator\n")
    log.write(f"CPU Count: {os.cpu_count()}\n")
    rmsd_distance_matrix = rdkitutils.rmsd_matrix_parallel(ens_mols)

    log.write(f"Exporting .csv distance matrix...\n")
    np.savetxt(f'{basename}.csv', rmsd_distance_matrix, delimiter=",")
    log.write(f"Done\n")

    to_delete = geometry.get_duplicates_rmsd_matrix(rmsd_distance_matrix, threshold=args.threshold)
    

    _coordinates = [ens_coordinates[i] for i in range(numconfs) if i not in to_delete]
    _elements = [ens_elements[i] for i in range(numconfs) if i not in to_delete]
    _header = [ens_header[i] for i in range(numconfs) if i not in to_delete]

    log.write(f"Conformers after deduplication: {len(_coordinates)}\n")
    log.write(f"{numconfs - len(_coordinates)} conformers were deleted\n")
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
        