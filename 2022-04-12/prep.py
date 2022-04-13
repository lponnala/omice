
import pandas as pd
from scipy.stats import ttest_ind

# X = pd.read_excel("/home/lponnala/repos/omics/2018-02-17/fort-test-Antonset.xlsx")
# X = X.iloc[:28,1:9]
# X = X.rename(columns={'simple functional name (July 2017)':'function'})

X = pd.read_excel("/home/lponnala/repos/omics/")

S = []
for _,x in X.iterrows():
    res = ttest_ind(x.iloc[2:5], x.iloc[5:8], equal_var=False, alternative='two-sided')
    S.append((x['function'], x['location'], res.statistic, res.pvalue))

S = pd.DataFrame(S, columns=['function','location','t-stat','p-value'])
S['infer'] = (S['p-value'] < 0.05).map({True: 'NOTEQ', False: 'EQUAL'})
print(S)
