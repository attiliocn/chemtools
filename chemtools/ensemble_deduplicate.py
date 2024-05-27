#!/usr/bin/env python3

import argparse
import time

from modules.xyzutils import read_xyz_ensemble, build_xyz_file
from modules import geometry

start_time_global = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+', help='XYZ Ensemble Files')
parser.add_argument('--threshold', type=float, default=.25, help='RMSD threshold for duplicate detection. Default is 0.25')
parser.add_argument('--parallel', action='store_true', help='Use the parallel version of the RMSD matrix calculator')
args = parser.parse_args()

log = open('deduplicate.log', mode='w', buffering=1)
for file in args.files:
    basename, extension = file.rsplit('.', 1)
    basename_updated = f"{basename}_deduplicated.{extension}"

    ensemble = read_xyz_ensemble(file)
    numconfs = len(ensemble)
    coordinates_all = [_['coordinates'] for _ in ensemble.values()]
    elements_all = [_['elements'] for _ in ensemble.values()]
    header_all = [_['header'] for _ in ensemble.values()]

    log.write(f"Current file: {basename}\n")
    log.write(f"Number of conformers: {numconfs}\n")
    start_time = time.time()

    if args.parallel:
        import os
        log.write("Using parallel RMSD calculator\n")
        log.write(f"CPU Count: {os.cpu_count()}\n")
        rmsd_distance_matrix = geometry.rmsd_matrix_parallel(coordinates_all)
    else:
        rmsd_distance_matrix = geometry.rmsd_matrix(coordinates_all)
    
    to_delete = geometry.get_duplicates_rmsd_matrix(rmsd_distance_matrix, threshold=args.threshold)

    _coordinates = [coordinates_all[i] for i in range(numconfs) if i not in to_delete]
    _elements = [elements_all[i] for i in range(numconfs) if i not in to_delete]
    _header = [header_all[i] for i in range(numconfs) if i not in to_delete]

    log.write(f"Conformers after deduplication: {len(_coordinates)}\n")
    log.write(f"{numconfs - len(_coordinates)} conformers were deleted\n")
    log.write(f"{len(_coordinates)} conformers will be written to {basename_updated}\n")
    end_time = time.time()
    elapsed_time = end_time - start_time
    log.write(f"Wall clock time: {elapsed_time:.2f} seconds\n\n")

    with open(basename_updated, mode='w') as f:
        for header, elements, coordinates in zip(_header, _elements, _coordinates):  
            _ = build_xyz_file(
                elements, 
                coordinates,
                header=header
            )
            f.write(_)

end_time_global = time.time()
elapsed_time_global = end_time_global - start_time_global
log.write(f"Total wall clock time: {elapsed_time_global:.2f} seconds\n")

log.close()
        