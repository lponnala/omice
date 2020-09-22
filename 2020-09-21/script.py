
import pandas as pd

A = pd.read_excel("F:/download/ATHENA-PGs-ABC1Ks-forLalit.xlsx",sheet_name=None,skiprows=1)

df = A['PGs'].set_index(['our interests', 'group', 'Accession', 'lab annotation'])
# df = A['17-ABC1Ks'].set_index(['our interests', 'group', 'Accession', 'lab annotation'])

if any(df.isna().sum(axis=0) == df.shape[0]):
    print(f"{sum(df.isna().sum(axis=0) == df.shape[0])} all-missing column(s) found!")
    df = df.loc[:,df.isna().sum(axis=0) < df.shape[0]]

if any(df.isna().sum(axis=1) == df.shape[1]):
    print(f"{sum(df.isna().sum(axis=1) == df.shape[1])} all-missing row(s) found!")
    df = df.loc[df.isna().sum(axis=1) < df.shape[1],:]

# fill missing values with zero
df = df.fillna(0)

# divide each row by its max (i.e. by the max across all tissues)
df_s = df.apply(lambda x: x/x.max(), axis=1)

# # min-max scaling of each row
# df_mm = df.apply(lambda x: (x-x.min())/(x.max()-x.min()), axis=1)

# standardize each row
df_z = df.apply(lambda x: (x-x.mean())/x.std(ddof=1), axis=1)


