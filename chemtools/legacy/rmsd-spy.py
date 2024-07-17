#!/usr/bin/env python3

import morfeus.utils
import spyrmsd.rmsd
from modules import xyzutils
import spyrmsd
import morfeus

ensemble = xyzutils.read_xyz_ensemble('ensemble.xyz')
numconfs = len(ensemble)
ens_elements = [_['elements'] for _ in ensemble.values()]
ens_atoms = [_['atomic_numbers'] for _ in ensemble.values()]
ens_coordinates = [_['coordinates'] for _ in ensemble.values()]
ens_header = [_['header'] for _ in ensemble.values()]
ens_cm = [
    morfeus.utils.get_connectivity_matrix(coordinates, elements) for elements,coordinates in zip(ens_elements,ens_coordinates)
]

rmsd = spyrmsd.rmsd.symmrmsd(
    coordsref=ens_coordinates[2],
    amref=ens_cm[2],
    apropsref=ens_atoms[2],
    coords=ens_coordinates[3],
    am=ens_cm[3],
    aprops=ens_atoms[3]
)

print(rmsd)