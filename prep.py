
import pandas as pd

X = pd.read_excel("/home/lponnala/repos/omics/2018-02-17/fort-test-Antonset.xlsx")
X = X.iloc[:,1:9]
X = X.rename(columns={'simple functional name (July 2017)':'function'})

