
import pandas as pd

data_file = "forLalit-3sets-forclusterproteome.xlsx"

# ~~ Set 1 ~~
X = pd.read_excel(data_file, sheet_name="set1-CLP", skiprows=1)

if X.shape[0] != len(X['Accession'].drop_duplicates()):
    raise Exception("duplicate accession numbers found!")

X = X.drop(['our interests','group','lab annotation'],axis=1)
print(X.isna().sum())
X = X.fillna(0).set_index('Accession').rename_axis(columns='Tissue')

# -- by protein --
P = X
# remove proteins that have all zero values (since they create problems for the correlation-based distance matrix)
P = P[(P == 0).sum(axis=1) != P.shape[1]]
# normalize each row by its max (i.e. by max across all tissues)
P = P.apply(lambda x: x/x.max(), axis=1)
# ensure no missing values
if P.isna().any().any():
    raise Exception("NaN values in P")
# write to file
P.reset_index().to_csv("set1_data_byProtein.csv",index=False)

# -- by tissue --
T = X.T
# remove tissues that have all zero values (since they create problems for the correlation-based distance matrix)
T = T[(T == 0).sum(axis=1) != T.shape[1]]
# normalize each row by its max (i.e. by max across all proteins)
T = T.apply(lambda x: x/x.max(), axis=1)
# ensure no missing values
if T.isna().any().any():
    raise Exception("NaN values in T")
# write to file
T.reset_index().to_csv("set1_data_byTissue.csv",index=False)


# ~~ Set 2 ~~
X = pd.read_excel(data_file, sheet_name="set2-PG-ABC", skiprows=1)

if X.shape[0] != len(X['Accession'].drop_duplicates()):
    raise Exception("duplicate accession numbers found!")

X = X.drop(['our interests','group','lab annotation'],axis=1)
print(X.isna().sum())
X = X.fillna(0).set_index('Accession').rename_axis(columns='Tissue')

# from scipy.cluster.hierarchy import linkage
# linked = linkage(X, method='average', metric='correlation')

# -- by protein --
P = X
# remove proteins that have all zero values (since they create problems for the correlation-based distance matrix)
P = P[(P == 0).sum(axis=1) != P.shape[1]]
# normalize each row by its max (i.e. by max across all tissues)
P = P.apply(lambda x: x/x.max(), axis=1)
# ensure no missing values
if P.isna().any().any():
    raise Exception("NaN values in P")
# write to file
P.reset_index().to_csv("set2_data_byProtein.csv",index=False)

# -- by tissue --
T = X.T
# remove tissues that have all zero values (since they create problems for the correlation-based distance matrix)
T = T[(T == 0).sum(axis=1) != T.shape[1]]
# normalize each row by its max (i.e. by max across all proteins)
T = T.apply(lambda x: x/x.max(), axis=1)
# ensure no missing values
if T.isna().any().any():
    raise Exception("NaN values in T")
# write to file
T.reset_index().to_csv("set2_data_byTissue.csv",index=False)


# ~~ Set 3 ~~
X = pd.read_excel(data_file, sheet_name="set3-TRAP", skiprows=1)

if X.shape[0] != len(X['Accession'].drop_duplicates()):
    raise Exception("duplicate accession numbers found!")

X = X.drop(['our interests','group','lab annotation'],axis=1)
print(X.isna().sum())
X = X.fillna(0).set_index('Accession').rename_axis(columns='Tissue')

P = X / X.max(axis=0)
if any(X.max(axis=1) == 0):
    T = X[X.max(axis=1) > 0].T / X[X.max(axis=1) > 0].max(axis=1)

if P.isna().any().any():
    raise Exception("NaN values in P")

if T.isna().any().any():
    raise Exception("NaN values in T")

P.reset_index().to_csv("set3_data_byProtein.csv",index=False)
T.reset_index().to_csv("set3_data_byTissue.csv",index=False)

