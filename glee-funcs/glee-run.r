source("glee-r-funcs.r")

# -- user-defined --

# XLFILE = "../glee/short.Cooper_147_vs_689.xls"
XLFILE = "Cooper_147_vs_689.xls"
# XLFILE = "wt-ClpS.xlsx"

# number of replicates
nA = 3
nB = 3

# specify rowIndex and/or colIndex if errors occur
Data = read.xlsx(XLFILE, sheetIndex=1, rowIndex=NULL, colIndex=NULL, as.data.frame=TRUE, header=TRUE)

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
