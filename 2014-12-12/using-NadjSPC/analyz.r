
# # Compare with the previous (2014 Feb) dataset

# library(xlsx)

# FILE = "data.xlsx"

# D = read.xlsx(FILE, sheetIndex=1, rowIndex=1:1881, colIndex=1:7, header=T)
# Dold = read.xlsx("../2014-02/raw-data.xls", sheetIndex=1, rowIndex=1:1881, colIndex=1:7, header=T)

# save(D, file="D.RData")
# save(Dold, file="Dold.RData")


# ==== Differential expression ====

source("glee-r-funcs.r")

# -- user-defined --

# library(xlsx)
# D = read.xlsx("data.xlsx", sheetIndex=1, rowIndex=1:1882, colIndex=1:7, header=T)
# save(D, file="D.RData")

load("D.RData")
Data = D

# number of replicates
nA = 3
nB = 3

# options
fit_type = "cubic"
num_iter = 10000
num_digits = 4
outfile = "diff-exp.xlsx"
# ------------------

if (!data_ok(Data,nA,nB)) {
	msg = paste("check input data! spreadsheet must have",
			"(1) the right number of columns",
			"(2) positive finite values",
			"\n",sep="\n")
	stop(msg)
}

Prot = as.character(Data[,1])
A = as.matrix(Data[,1+(1:nA)])
B = as.matrix(Data[,1+nA+(1:nB)])

m = fit_model(A, B, fit_type)
model_fit_plots(m, outfile="fitplots.png")
stn_pval = calc_stn_pval(A, B, m, num_iter)
stn_pval_plots(stn_pval, outfile="stn-pval.png")
diff_exp_table(stn_pval, Prot, num_digits, outfile)
