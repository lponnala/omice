
# ==== GLEE ====

source("glee-r-funcs.r")

# -- user-defined --

# library(xlsx)
# D = read.xlsx("data-adjSPC.xlsx", sheetIndex=1, rowIndex=1:1882, colIndex=1:7, header=T, stringsAsFactors=T)
# save(D, file="D.RData")

load("D.RData")
for (j in 2:7) {
	ind=which(is.na(D[,j]))
	D[ind,j]=0
}
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


# ==== QSPEC ====

# create the data matrix
load("D.RData")
for (j in 2:7) {
	ind=which(is.na(D[,j]))
	D[ind,j]=0
}
for (j in 1:dim(D)[2]) {
	D[,j] = as.character(D[,j])
}
# insert a column of ones as a dummy for length
D =cbind(D[,1], as.character(rep(1,dim(D)[1])), D[,-1], stringsAsFactors=F)
outfile = "datamatrix.txt"
cat("Protein\tLength\t0\t0\t0\t1\t1\t1\n", file=outfile, sep="", append=F)
for (i in 1:dim(D)[1]) {
	cat(unlist(D[i,]), file=outfile, sep="\t", append=T)
	cat("\n", file=outfile, append=T)
}

# // LINUX VERSION //
# run these commands in the Linux terminal window
# (working dir: /home/lalit/Code/R/omics/2014-12-12/using-adjSPC)
# $ LD_LIBRARY_PATH=/usr/local/lib
# $ export LD_LIBRARY_PATH
# $ ../../qprot-linux/qprot_1.3.0/qspec-param datamatrix.txt 5000 20000 0
# $ ../../qprot-linux/qprot_1.3.0/getfdr datamatrix.txt_qspec

# // WIN64 VERSION //
# run these commands in the Cygwin terminal window
# (working dir: ~/CODE/R/omics/2014-12-12/using-adjSPC)
# ../../qprot-win64/qprot_1.3.0/qspec-param.exe datamatrix.txt 5000 20000 0
# ../../qprot-win64/qprot_1.3.0/getfdr.exe datamatrix.txt_qspec

# NOTE:
# I have checked both versions above and they yield identical output, so either could be used
# The win64 and win32 versions yield identical output when tested with dummy data

# clean up the output (remove the "dummy" length column and sort in FDR order)
R = read.table("datamatrix.txt_qspec_fdr", header=T, sep="\t")
R = R[,-2]
names(R)[2:7] = c( paste("wt-",1:3,sep=""), paste("t1t2-",1:3,sep="") )
R = R[order(R$fdr),]
library(xlsx)
write.xlsx(R, file="qspec_results.xlsx", col.names=T, row.names=F)


# ==== Compare GLEE and QSPEC results ====

GLEE_results_file = "diff-exp.xlsx"
GLEE_pVal_cutoff = 0.01
QSPEC_results_file = "qspec_results.xlsx"
QSPEC_FDR_cutoff = 0.1

G = read.xlsx(GLEE_results_file, 1, header=T)
Q = read.xlsx(QSPEC_results_file, 1, header=T)

G_signif_ind = which(as.numeric(as.character(G$pVal)) < GLEE_pVal_cutoff)
G_signif_prot = as.character(G[G_signif_ind,1])
cat("\n", "# signif proteins as per GLEE = ", length(G_signif_prot), "\n")

Q_signif_ind = which(Q$fdr < QSPEC_FDR_cutoff)
Q_signif_prot = as.character(Q[Q_signif_ind,1])
cat("\n", "# signif proteins as per QSPEC = ", length(Q_signif_prot), "\n")

cat("\n", "number of signif proteins identified by either = ", length(union(G_signif_prot, Q_signif_prot)))
cat("\n", "number of signif proteins identified by both = ", length(intersect(G_signif_prot, Q_signif_prot)))
