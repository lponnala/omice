
import pandas as pd

data_file = "forLalit-3sets-forclusterproteome.xlsx"
sheets = ['set1-CLP','set2-PG-ABC','set3-TRAP']

for i,sheet in enumerate(sheets):
    print("\n",f"~~ processing {sheet} ~~","\n",sep="")
    X = pd.read_excel(data_file, sheet_name=sheet, skiprows=1)

    if X.shape[0] != len(X['Accession'].drop_duplicates()):
        raise Exception("duplicate accession numbers found!")

    X = X.drop(['our interests','group','lab annotation'],axis=1)
    # print(X.isna().sum())
    X = X.fillna(0).set_index('Accession').rename_axis(columns='Tissue')

    # -- by protein --
    P = X
    # remove proteins that have all zero values (since they create problems for the correlation-based distance matrix)
    if any((P == 0).sum(axis=1) == P.shape[1]):
        print(f"dropping {sum((P == 0).sum(axis=1) == P.shape[1])} proteins")
        P = P[(P == 0).sum(axis=1) != P.shape[1]]
    # normalize each row by its max (i.e. by the max across all tissues)
    P = P.apply(lambda x: x/x.max(), axis=1)
    # ensure no missing values
    if P.isna().any().any():
        raise Exception("missing values in P")
    # write to file
    P.reset_index().to_csv(f"set{i+1}_data_byProtein.csv",index=False)

    # -- by tissue --
    T = X.T
    # remove tissues that have all zero values (since they create problems for the correlation-based distance matrix)
    if any((T == 0).sum(axis=1) == T.shape[1]):
        print(f"dropping {sum((T == 0).sum(axis=1) == T.shape[1])} tissues")
        T = T[(T == 0).sum(axis=1) != T.shape[1]]
    # normalize each row by its max (i.e. by the max across all proteins)
    T = T.apply(lambda x: x/x.max(), axis=1)
    # ensure no missing values
    if T.isna().any().any():
        raise Exception("missing values in T")
    # write to file
    T.reset_index().to_csv(f"set{i+1}_data_byTissue.csv",index=False)

