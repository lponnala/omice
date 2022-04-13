
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats import multicomp, multitest

D = pd.read_excel("/home/lponnala/repos/omics/2022-04-12/forlalit-chloroplast-wt-prep1prep2.xlsx")
X, ofile = D.iloc[:,[0,4,5,6,7,8,9]], 'all.xlsx'
X, ofile = D[(D.iloc[:,[10,11,12]].sum(axis=1) > 6) | (D.iloc[:,[13,14,15]].sum(axis=1) > 6)].iloc[:,[0,4,5,6,7,8,9]], 'subset.xlsx'

X.shape
X.dtypes

S = []
for _,x in X.iterrows():
    res = ttest_ind(x.iloc[1:4], x.iloc[4:7], equal_var=False, alternative='two-sided')
    S.append((x['protein'], res.statistic, res.pvalue))

S = pd.DataFrame(S, columns=['protein','t-stat','p-value'])
S['infer'] = (S['p-value'] < 0.05).map({True: 'NOTEQ', False: 'EQUAL'})
S['infer'].value_counts()
print(S)

y = multitest.multipletests(S['p-value'], alpha=0.1, method='bonferroni')
y[1]

multicomp.tukeyhsd()

# # -- verify previous result --
# # X = pd.read_excel("/home/lponnala/repos/omics/2018-02-17/fort-test-Antonset.xlsx")
# # X = X.iloc[:28,1:9]
# # X = X.rename(columns={'simple functional name (July 2017)':'function'})
# S = []
# for _,x in X.iterrows():
#     res = ttest_ind(x.iloc[2:5], x.iloc[5:8], equal_var=False, alternative='two-sided')
#     S.append((x['function'], x['location'], res.statistic, res.pvalue))
# S = pd.DataFrame(S, columns=['function','location','t-stat','p-value'])
# S['infer'] = (S['p-value'] < 0.05).map({True: 'NOTEQ', False: 'EQUAL'})
# print(S)
