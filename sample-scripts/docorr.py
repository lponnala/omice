from __future__ import division
import xlrd
import numpy as np
from scipy.stats import pearsonr
import itertools

FILE = 'correlation-wt-clps-clpc.xls'
wb = xlrd.open_workbook(FILE)
sh = wb.sheet_by_index(0)

def print_corr(D):
	for i,j in itertools.combinations(range(len(D)),2):
		I = (D[i]+D[j])>0
		(corr,pvalue) = pearsonr(D[i][I],D[j][I])
		print 'rep%d,rep%d\t%d\t%g\t%g' % (i+1, j+1, sum(I), corr, pvalue)

P = sh.col_values(0,1)
WT=[]
for j in range(1,1+3):
	WT.append(sh.col_values(j,1))
# WT = np.column_stack(WT)
WT = np.asarray(WT)
print 'WT\tNUM_PROTEINS\tCORRELATION\tPVALUE'
print_corr(WT)

clpS=[]
for j in range(4,4+3):
	clpS.append(sh.col_values(j,1))
# clpS = np.column_stack(clpS)
clpS = np.asarray(clpS)
print 'clpS\tNUM_PROTEINS\tCORRELATION\tPVALUE'
print_corr(clpS)

clpC=[]
for j in range(7,7+3):
	clpC.append(sh.col_values(j,1))
# clpC = np.column_stack(clpC)
clpC = np.asarray(clpC)
print 'clpC\tNUM_PROTEINS\tCORRELATION\tPVALUE'
print_corr(clpC)

