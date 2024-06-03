#!/usr/bin/env python3

import argparse
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score

import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument('distance_matrix', nargs='+', help='.csv distance matrix')
parser.add_argument('--method', default='distance_threshold', choices=['distance_threshold','n_confs'], help='Method to generate the clusters')
parser.add_argument('--num-clusters', type=int, default=None, help='Do not test for the best number of clusters. Use INT number of clusters instead')
parser.add_argument('--max-clusters', type=int, default=50, help='Maximum number of clusters to test (--method=n_confs only)')
parser.add_argument('--distance-threshold', type=float, default=1.0, help=' Maximum distance inside each cluster in A (--method=distance_threshold only)')
parser.add_argument('--triangular', action='store_true', help='The RMSD matrix is in the (lower) triangular form')
args = parser.parse_args()

for file in args.distance_matrix:
    basename = file.rsplit('.',1)[0]

    dist_matrix =  pd.read_csv(file, index_col=None, header=None)
    
    if args.triangular:
        dist_matrix += dist_matrix.T # convert matrix from lower triangular to symmetric

    if args.method == 'n_confs':
        if args.num_clusters:
            n_clusters = args.num_clusters
            clustering = AgglomerativeClustering(
                n_clusters=n_clusters,
                compute_full_tree=True,
                distance_threshold=None,
                metric='precomputed',
                linkage='average'
            )
        else:
            n_clusters = args.max_clusters
            method_score = {}
            for n_clusters in range(2,n_clusters):
                clustering = AgglomerativeClustering(
                    n_clusters=n_clusters,
                    compute_full_tree=True,
                    distance_threshold=None,
                    metric='precomputed',
                    linkage='average'
                )
                predictions = clustering.fit_predict(dist_matrix)
                score = {
                    'silhouette': silhouette_score(dist_matrix, predictions), # higher
                    'calinski': calinski_harabasz_score(dist_matrix, predictions), # higher
                    'davies': davies_bouldin_score(dist_matrix, predictions) # lower
                }
                method_score[n_clusters] = score
            
            method_score_df = pd.DataFrame(method_score).T

            fig, [ax1,ax2,ax3] = plt.subplots(nrows=3)
            sns.lineplot(ax=ax1, x=method_score_df.index, y=method_score_df['silhouette'], label='silhouette')
            sns.lineplot(ax=ax2, x=method_score_df.index, y=method_score_df['calinski'], color='orange', label='calinsky')
            sns.lineplot(ax=ax3, x=method_score_df.index, y=method_score_df['davies'], color='red', label='davies')
            for ax in [ax1, ax2, ax3]:
                ax.set_xlabel('')
                ax.set_ylabel('')
            fig.tight_layout()
            plt.savefig('cluster-metrics.svg')
            n_clusters = int(method_score_df['calinski'].idxmax())
    
        clustering = AgglomerativeClustering(
                n_clusters=n_clusters,
                compute_full_tree=True,
                distance_threshold=None,
                metric='precomputed',
                linkage='average'
            )
        
    elif args.method == 'distance_threshold':
        clustering = AgglomerativeClustering(
            n_clusters=None,
            compute_full_tree=True,
            distance_threshold=args.distance_threshold,
            metric='precomputed',
            linkage='average'
        )
    
    predictions = clustering.fit_predict(dist_matrix)
    n_clusters = clustering.n_clusters_
    print(f"Total number of clusters: {n_clusters}")

    with open(f'{basename}.sh', mode='w') as f:
        f.write(f'mkdir {basename}\n')
        f.write(f'mv {basename}* {basename}\n')
        f.write(f'mv cluster-metrics.svg {basename}\n' )
        f.write(f'cd {basename}\n')
        f.write(f'split_conformer_ensemble.py {basename}.xyz &> /dev/null\n')
        f.write(f'mkdir cluster-{{1..{n_clusters}}}\n')
        for i, j in enumerate(predictions):
            f.write(f'mv *_{i+1}.xyz cluster-{j+1}/\n')
        f.write('mkdir ensemble && for f in cluster-*/; do find $f -type f -name "*.xyz" | sort -V | head -n 1; done | while read -r f; do cp $f ensemble; done\n')
        f.write('cd ..\n')

    
