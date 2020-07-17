
import pandas as pd

data_file = "forLalit-3sets-forclusterproteome.xlsx"

X = pd.read_excel(data_file, sheet_name="set1-CLP", skiprows=1)

if X.shape[0] != len(X['Accession'].drop_duplicates()):
    raise Exception("duplicate accession numbers found!")

X = X.drop(['our interests','group','lab annotation'],axis=1)
print(X.isna().sum())
X = X.fillna(0)



