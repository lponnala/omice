
# https://www.analyticsvidhya.com/blog/2019/05/beginners-guide-hierarchical-clustering/
# https://scikit-learn.org/stable/modules/clustering.html#hierarchical-clustering

# ----------------------------------------
# // Heirarchical Clustering //
# see 2014-11-24/analyz.R

# ~~ prepare the data ~~
import pandas as pd
clust_datafile = "ABC1K_proteomicsDB_expression-clustering-forLalit.xlsx"

D = pd.read_excel(clust_datafile)
# len(D['Accession'].unique())
# len(D['Tissue (27 types)'].dropna().values)
# D['Tissue'].value_counts()
D = D.iloc[:,:4]

# per Klaas' email, normalize each protein by the max expression value across 27 tissue types
maxValue = D.groupby('Accession')['Average Normalized Intensity'].max().reset_index().rename({'Average Normalized Intensity':'Max Intensity'},axis=1)
P = pd.merge(D, maxValue, how='left', on='Accession').assign(Intensity = lambda X: X['Average Normalized Intensity']/X['Max Intensity'])
P = P.pivot(index='Accession', columns='Tissue', values='Intensity').fillna(0)

# per Klaas' second email, normalize each tissue by the max expression value across all proteins
maxValue = D.groupby('Tissue')['Average Normalized Intensity'].max().reset_index().rename({'Average Normalized Intensity':'Max Intensity'},axis=1)
T = pd.merge(D, maxValue, how='left', on='Tissue').assign(Intensity = lambda X: X['Average Normalized Intensity']/X['Max Intensity'])
T = T.pivot(index='Tissue', columns='Accession', values='Intensity').fillna(0)

# save the normalized-to-max data to be used for clustering
P.reset_index().to_csv("cluster_data_byProtein.csv",index=False)
T.reset_index().to_csv("cluster_data_byTissue.csv",index=False)

# # ~~ draw dendrogram to identify the number of clusters ~~
# from scipy.cluster.hierarchy import dendrogram, linkage
# from matplotlib import pyplot as plt
# from sklearn.cluster import AgglomerativeClustering
# # !! can't mark rectanges around the clusters, so go with the R version in Rscript.R !!
# linked = linkage(P, method='average', metric='correlation')
# labelList = P.index.values
# plt.figure(figsize=(10, 7))
# dendrogram(linked,
#             # color_threshold=0.8,
#             truncate_mode='lastp', p = 17,
#             orientation='top',
#             labels=labelList,
#             distance_sort='descending',
#             show_leaf_counts=True, leaf_rotation=90)
# plt.tight_layout()
# plt.show()
# # ~~ find the contents of each cluster (matches with the R version in Rscript.R) ~~
# num_clust = 4
# cluster = AgglomerativeClustering(n_clusters=num_clust, affinity='correlation', linkage='average')
# y = cluster.fit_predict(P)
# for k in range(num_clust):
#     print(f"cluster {k+1}")
#     print(P[y == k].index)


# ----------------------------------------
# // Differential Expression //
# see 2018-02-17/compare-across-funcs.R

import pandas as pd

diffexp_datafile = "pgs-elena-forlalit-Glee_annovatest.xlsx"

D = pd.read_excel(diffexp_datafile)
D = D.iloc[:,:7].fillna(0)
D.columns = D.columns.map(lambda x: x.lower().replace(" pg replicate ","_")).map(lambda x: x.replace(" ","_"))

# prepare data for glee
D.loc[:,['protein_id','wt_1','wt_2','k1_1','k1_2']].to_csv("data_glee_wt-k1.csv", index=False)
D.loc[:,['protein_id','k1_1','k1_2','k6_1','k6_2']].to_csv("data_glee_k1-k6.csv", index=False)
D.loc[:,['protein_id','wt_1','wt_2','k6_1','k6_2']].to_csv("data_glee_wt-k6.csv", index=False)

# prepare data for anova
D = D.melt(id_vars='protein_id', var_name='gen_rep', value_name='val')
D['gen'] = D['gen_rep'].map(lambda x: x.split('_')[0])
D['rep'] = D['gen_rep'].map(lambda x: 'rep' + x.split('_')[1])
D.drop('gen_rep',axis=1,inplace=True)
D.groupby(['gen','rep']).size()
D.to_csv("diffexp_anova_data.csv", index=False)
