
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats import weightstats, multitest, multicomp

D = pd.read_excel("/home/lponnala/repos/omics/2022-04-12/forlalit-chloroplast-wt-prep1prep2.xlsx")
ttest_from = ['scipy','statsmodels'][1]
dataset = ['all','subset'][0]

if dataset == 'all':
    X = D.iloc[:,[0,4,5,6,7,8,9]]
elif dataset == 'subset':
    X = D[(D.iloc[:,[10,11,12]].sum(axis=1) > 6) | (D.iloc[:,[13,14,15]].sum(axis=1) > 6)].iloc[:,[0,4,5,6,7,8,9]]

X.shape
X.dtypes

S = []
for _,x in X.iterrows():
    if ttest_from == 'scipy':
        res = ttest_ind(x.iloc[1:4], x.iloc[4:7], equal_var=False, alternative='two-sided')
        S.append((x['protein'], res.statistic, res.pvalue))
    elif ttest_from == 'statsmodels':
        res = weightstats.ttest_ind(x.iloc[1:4], x.iloc[4:7], usevar='unequal', alternative='two-sided')
        S.append((x['protein'], res[0], res[1]))

S = pd.DataFrame(S, columns=['protein','t-stat','p-value'])
S['inference'] = (S['p-value'] < 0.05).map({True: 'NOTEQ', False: 'EQUAL'})
S['inference'].value_counts()
# print(S)

for method in ['bonferroni','sidak','holm-sidak','holm','simes-hochberg','hommel','fdr_bh','fdr_by','fdr_tsbh','fdr_tsbky'][:1]:
    print("\n",f".. {method} ..", sep="")
    y = multitest.multipletests(S['p-value'], alpha=0.1, method=method)
    S['corrected-p-value'] = y[1]
    S['corrected-inference'] = (S['corrected-p-value'] < 0.05).map({True: 'NOTEQ', False: 'EQUAL'})
    S['corrected-inference'].value_counts()

# multicomp.tukeyhsd()

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
