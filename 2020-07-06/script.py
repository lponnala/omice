
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
Y = pd.merge(D, maxValue, how='left', on='Accession').assign(Intensity = lambda X: X['Average Normalized Intensity']/X['Max Intensity'])
Y = Y.pivot(index='Accession', columns='Tissue', values='Intensity').fillna(0)


# # ----------------------------------------
# # // Differential Expression //
# 
# diffexp_datafile = "pgs-elena-forlalit-Glee_annovatest.xlsx"
# D = pd.read_excel(diffexp_datafile)

