import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from modules.xyzutils import read_xyz_ensemble, build_xyz_file
from modules.geometry import rmsd_matrix_compare

parser = argparse.ArgumentParser()
parser.add_argument('ensemble1', help='XYZ Ensemble no 1')
parser.add_argument('ensemble2', help='XYZ Ensemble no 2')
parser.add_argument('--no-align', action='store_false', help='Ignore matrix aligment')
args = parser.parse_args()

ensemble_a = read_xyz_ensemble(args.ensemble1)
coordinates_a = [_['coordinates'] for _ in ensemble_a.values()]
ensemble_b = read_xyz_ensemble(args.ensemble2)
coordinates_b = [_['coordinates'] for _ in ensemble_b.values()]

rmsd_matrix = rmsd_matrix_compare(coordinates_a, coordinates_b, align=args.no_align)
rmsd_matrix = rmsd_matrix.round(3)

np.savetxt("rmsd-compare.csv", rmsd_matrix, delimiter=",")

plt.close('all')
fig, ax = plt.subplots(dpi=300)
heatmap = sns.heatmap(
    data=rmsd_matrix
)
fig.tight_layout()
plt.savefig('rmsd-compare-rmsd.png')

plt.close('all')
fig, ax = plt.subplots(dpi=300)
cluster = sns.clustermap(
    data=rmsd_matrix
)
fig.tight_layout()
plt.savefig('rmsd-compare-cluster.png')