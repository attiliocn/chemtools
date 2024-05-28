#! /usr/bin/env python3

import argparse
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdMolAlign

import matplotlib.pyplot as plt
import seaborn as sns

from modules.xyzutils import read_xyz_ensemble, read_xyz_file
from modules.rdkitutils import convert_coordinates_to_mols

parser = argparse.ArgumentParser()
parser.add_argument('ensemble1', help='XYZ Ensemble no 1')
parser.add_argument('ensemble2', help='XYZ Ensemble no 2')
args = parser.parse_args()

log = open('rmsd-compare.log', mode='w', buffering=1)
log.write('Conformational Ensemble Comparison starting...\n')
log.write('Reading ensembles\n')

ensemble_a = read_xyz_ensemble(args.ensemble1)
elements_a = [_['elements'] for _ in ensemble_a.values()]
coordinates_a = [_['coordinates'] for _ in ensemble_a.values()]
mols_a = [convert_coordinates_to_mols(ele, coords) for ele, coords in zip(elements_a, coordinates_a)]
mols_a = np.array(mols_a)

ensemble_b = read_xyz_ensemble(args.ensemble2)
elements_b = [_['elements'] for _ in ensemble_b.values()]
coordinates_b = [_['coordinates'] for _ in ensemble_b.values()]
mols_b = [convert_coordinates_to_mols(ele, coords) for ele, coords in zip(elements_b, coordinates_b)]
mols_b = np.array(mols_b)

log.write('Finished importing the ensembles\n')

log.write(f"Ensemble A has {len(mols_a)} molecules\n")
log.write(f"Ensemble B has {len(mols_b)} molecules\n")
log.write(f"A ({len(mols_a)},{len(mols_b)}) matrix of RMSD will be calculated\n")

all_mols = np.concatenate([mols_a, mols_b])

log.write('Determining MCS (Maximum Common Substructure)\n')

mcs = rdFMCS.FindMCS(all_mols)
pattern = Chem.MolFromSmarts(mcs.smartsString)

log.write('MCS was determined using all atoms\n')
log.write(f"MCS has {mcs.numAtoms} atoms\n")
log.write(f"MCS SMARTS Pattern is {mcs.smartsString}\n\n")
log.write('Removing hydrogens from the SMARTS patterns\n')
pattern = Chem.RemoveAllHs(pattern)
mcs_heavy = Chem.MolToSmarts(pattern)
log.write(f'MCS SMARTS Pattern (heavy-atoms) is {mcs_heavy}\n\n')

log.write('Calculating heavy-atom RMSD Matrix (it may take a while)\n')

rmsd_matrix = np.zeros((len(mols_a), len(mols_b)))
for i in range(rmsd_matrix.shape[0]):
    mol_ref = mols_a[i]
    ref_match = mol_ref.GetSubstructMatch(pattern)
    for j in range(rmsd_matrix.shape[1]):
        probe_mol = mols_b[j]
        mv = probe_mol.GetSubstructMatch(pattern)
        rmsd = rdMolAlign.AlignMol(probe_mol, mol_ref, atomMap=list(zip(mv,ref_match)))
        rmsd_matrix[i, j] = rmsd

with open(f"{args.ensemble2.rsplit('.')[0]}_aligned.xyz", mode='w') as f:
    for mol in mols_b:
        f.write(Chem.MolToXYZBlock(mol))

log.write('Finished\n')
log.write('\n\n')
log.close()

np.savetxt("rmsd-compare.csv", rmsd_matrix, delimiter=",")

plt.close('all')
fig, ax = plt.subplots(dpi=300)
heatmap = sns.heatmap(
    data=rmsd_matrix
)
heatmap.set_xlabel(f"{args.ensemble2}")
heatmap.set_ylabel(f"{args.ensemble1}")
fig.tight_layout()
plt.savefig('rmsd-compare-rmsd.png')

plt.close('all')
fig, ax = plt.subplots(dpi=300)
cluster = sns.clustermap(
    data=rmsd_matrix
)
fig.tight_layout()
plt.savefig('rmsd-compare-cluster.png')