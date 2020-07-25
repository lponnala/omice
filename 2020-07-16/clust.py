
# This is an alternative (fully in Python) to clust.R

# ~~ draw dendrogram ~~
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt

from itertools import product
import seaborn as sns

for s,t in product():
    print("\n",f"processing {s} {t}","\n",sep="")
    X = pd.read_csv(f"{s}_data_{t}.csv")
    X = X.set_index({'byProtein': 'Accession', 'byTissue': 'Tissue'}[t])
    sns.clustermap(X, method='average', metric='correlation', row_cluster=True, col_cluster=False)
    plt.show()


for s in ('set1','set2','set3'):
    for t in ('byProtein','byTissue'):
        print("\n",f"processing {s} {t}","\n",sep="")
        X = pd.read_csv(f"{s}_data_{t}.csv")
        X = X.set_index({'byProtein': 'Accession', 'byTissue': 'Tissue'}[t])
        linked = linkage(X, method='average', metric='correlation')
        labelList = X.index.values
        plt.figure(figsize=(10, 7))
        dendrogram(linked,
                    color_threshold=0,
                    above_threshold_color='blue',
                    orientation='top',
                    labels=labelList,
                    distance_sort='descending',
                    show_leaf_counts=True, leaf_rotation=90)
        plt.title(f"{s}: Hierarchical clustering {t}")
        plt.tight_layout()
        plt.savefig(f"pyfigs/{s}_clusters_{t}.png")
        # plt.show()

# # ~~ find the contents of each cluster ~~
# from sklearn.cluster import AgglomerativeClustering
# num_clust = 4
# cluster = AgglomerativeClustering(n_clusters=num_clust, affinity='correlation', linkage='average')
# y = cluster.fit_predict(P)
# for k in range(num_clust):
#     print(f"cluster {k+1}")
#     print(P[y == k].index)

# from scipy.cluster.hierarchy import cut_tree
# from collections import Counter
# ctree = cut_tree(linked, n_clusters=12)
# Counter(ctree.squeeze())
