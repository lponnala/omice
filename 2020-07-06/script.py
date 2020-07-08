
# https://www.analyticsvidhya.com/blog/2019/05/beginners-guide-hierarchical-clustering/
# https://scikit-learn.org/stable/modules/clustering.html#hierarchical-clustering

# ----------------------------------------
# // Heirarchical Clustering //
# see 2014-11-24/analyz.R

import pandas as pd
clust_datafile = "ABC1K_proteomicsDB_expression-clustering-forLalit.xlsx"
D = pd.read_excel(clust_datafile)

# len(D['Accession'].unique())
# len(D['Tissue (27 types)'].dropna().values)
# D['Tissue'].value_counts()

D = D.iloc[:,:4]
maxValue = D.groupby('Accession')['Average Normalized Intensity'].max().reset_index().rename({'Average Normalized Intensity':'Max Intensity'},axis=1)
X = pd.merge(D, maxValue, how='left', on='Accession').assign(Intensity = lambda X: X['Average Normalized Intensity']/X['Max Intensity'])
X = X.pivot(index='Accession', columns='Tissue', values='Intensity').fillna(0)
X.reset_index().to_csv("cluster_data.csv",index=False)

from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt

linked = linkage(X, method='average', metric='correlation')
labelList = X.index.values
plt.figure(figsize=(10, 7))
dendrogram(linked,
            orientation='top',
            labels=labelList,
            distance_sort='descending',
            show_leaf_counts=True)
plt.show()

# import numpy as np
# from scipy.cluster.hierarchy import dendrogram, linkage
# from matplotlib import pyplot as plt
# 
# X = np.array([[5,3],
#     [10,15],
#     [15,12],
#     [24,10],
#     [30,30],
#     [85,70],
#     [71,80],
#     [60,78],
#     [70,55],
#     [80,91],])
# linked = linkage(X, 'single')
# labelList = range(1, 11)
# plt.figure(figsize=(10, 7))
# dendrogram(linked,
#             orientation='top',
#             labels=labelList,
#             distance_sort='descending',
#             show_leaf_counts=True)
# plt.show()


# # ----------------------------------------
# # // Differential Expression //
# # see 2018-02-17/compare-across-funcs.R
# 
# diffexp_datafile = "pgs-elena-forlalit-Glee_annovatest.xlsx"
# D = pd.read_excel(diffexp_datafile)
# D = D.iloc[:,:7]
# D.columns = D.columns.map(lambda x: x.lower().replace(" pg replicate ","_")).map(lambda x: x.replace(" ","_"))

