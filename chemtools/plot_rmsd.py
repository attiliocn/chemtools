#! /usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-f','--matrix', help='RMSD matrix in csv format')
parser.add_argument('--triangular', action='store_true', help='The RMSD matrix is in the (lower) triangular form')
parser.add_argument('--linkage', default='average', help='Linkage Method for Hierarchical Clustering. See <https://seaborn.pydata.org/generated/seaborn.clustermap.html> for details')
args = parser.parse_args()

rmsd_matrix = pd.read_csv(args.matrix, header=None, index_col=None)

if args.triangular:
    rmsd_matrix += rmsd_matrix.T

plt.close('all')
fig, ax = plt.subplots(dpi=300)
heatmap = sns.heatmap(
    data=rmsd_matrix
)
heatmap.set_xlabel(f"")
heatmap.set_ylabel(f"")
fig.tight_layout()
plt.savefig('plot_rmsd.png')

plt.close('all')
fig, ax = plt.subplots(dpi=300)
cluster = sns.clustermap(
    data=rmsd_matrix,
    method=args.linkage
)
fig.tight_layout()
plt.savefig('plot_rmsd-cluster.png')