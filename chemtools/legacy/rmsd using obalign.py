#!/usr/bin/env python3
from openbabel import openbabel
from modules import xyzutils
from modules import babelutils

ensemble = xyzutils.read_xyz_ensemble('ensemble.xyz')
mols = babelutils.convert_ensemble_to_mols(ensemble)

obalign = openbabel.OBAlign(False,True)

obalign.SetRefMol(mols[2])
obalign.SetTargetMol(mols[3])

obalign.Align()
print(obalign.GetRMSD())
