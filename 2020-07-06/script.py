
# https://www.analyticsvidhya.com/blog/2019/05/beginners-guide-hierarchical-clustering/
# https://scikit-learn.org/stable/modules/clustering.html#hierarchical-clustering

# ----------------------------------------
# // Heirarchical Clustering //
# see 2014-11-24/analyz.R

import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
from sklearn.cluster import AgglomerativeClustering

clust_datafile = "ABC1K_proteomicsDB_expression-clustering-forLalit.xlsx"

D = pd.read_excel(clust_datafile)
# len(D['Accession'].unique())
# len(D['Tissue (27 types)'].dropna().values)
# D['Tissue'].value_counts()

# per Klaas' email, normalize each protein by the max expression value across 27 tissue types
D = D.iloc[:,:4]
maxValue = D.groupby('Accession')['Average Normalized Intensity'].max().reset_index().rename({'Average Normalized Intensity':'Max Intensity'},axis=1)
X = pd.merge(D, maxValue, how='left', on='Accession').assign(Intensity = lambda X: X['Average Normalized Intensity']/X['Max Intensity'])
X = X.pivot(index='Accession', columns='Tissue', values='Intensity').fillna(0)

# save the normalized-to-max data to be used for clustering
X.reset_index().to_csv("cluster_data.csv",index=False)

# draw dendrogram to identify the number of clusters
# !! can't mark rectanges around the clusters, so go with the R version in Rscript.R !!
linked = linkage(X, method='average', metric='correlation')
labelList = X.index.values
plt.figure(figsize=(10, 7))
dendrogram(linked,
            # color_threshold=0.8,
            truncate_mode='lastp', p = 17,
            orientation='top',
            labels=labelList,
            distance_sort='descending',
            show_leaf_counts=True, leaf_rotation=90)
plt.tight_layout()
plt.show()

# find the contents of each cluster (matches with the R version in Rscript.R)
num_clust = 4
cluster = AgglomerativeClustering(n_clusters=num_clust, affinity='correlation', linkage='average')
y = cluster.fit_predict(X)
for k in range(num_clust):
    print(f"cluster {k+1}")
    print(X[y == k].index)


# ----------------------------------------
# // Differential Expression //
# see 2018-02-17/compare-across-funcs.R

import pandas as pd

diffexp_datafile = "pgs-elena-forlalit-Glee_annovatest.xlsx"

D = pd.read_excel(diffexp_datafile)
D = D.iloc[:,:7].fillna(0)
D.columns = D.columns.map(lambda x: x.lower().replace(" pg replicate ","_")).map(lambda x: x.replace(" ","_"))

