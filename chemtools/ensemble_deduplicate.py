#!/usr/bin/env python3

import argparse

from modules.xyzutils import read_xyz_ensemble, build_xyz_file
from modules import geometry

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

    if args.parallel:
        log.write("Using parallel RMSD calculator\n")
        rmsd_distance_matrix = geometry.rmsd_matrix_parallel(coordinates_all)
    else:
        rmsd_distance_matrix = geometry.rmsd_matrix(coordinates_all)
    
    to_delete = geometry.get_duplicates_rmsd_matrix(rmsd_distance_matrix, threshold=args.threshold)

    _coordinates = [coordinates_all[i] for i in range(numconfs) if i not in to_delete]
    _elements = [elements_all[i] for i in range(numconfs) if i not in to_delete]
    _header = [header_all[i] for i in range(numconfs) if i not in to_delete]

    log.write(f"Conformers after deduplication: {len(_coordinates)}\n")
    log.write(f"{numconfs - len(_coordinates)} conformers were deleted\n")
    log.write(f"{len(_coordinates)} conformers will be written to {basename_updated}\n\n")

    with open(basename_updated, mode='w') as f:
        for header, elements, coordinates in zip(_header, _elements, _coordinates):  
            _ = build_xyz_file(
                elements, 
                coordinates,
                header=header
            )
            f.write(_)
log.close()
        