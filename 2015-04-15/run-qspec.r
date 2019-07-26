
# ==== QSPEC ====

# create the data matrix
load("A.RData")
for (j in 1:ncol(A)) {
	ind = which(is.na(A[,j]))
	A[ind,j] = 0
}

load("Acc.RData")

# comp = "prep-wt"; Data = cbind(Acc, A[,c(2,6,10,1,5,9)])
# comp = "oop-wt"; Data = cbind(Acc, A[,c(3,7,11,1,5,9)])
# comp = "triple-wt"; Data = cbind(Acc, A[,c(4,8,12,1,5,9)])
# comp = "prep-triple"; Data = cbind(Acc, A[,c(2,6,10,4,8,12)])
comp = "oop-triple"; Data = cbind(Acc, A[,c(3,7,11,4,8,12)])

for (j in 1:dim(Data)[2]) {
	Data[,j] = as.character(Data[,j])
}
# insert a column of ones as a dummy for length
D = cbind(Data[,1], as.character(rep(1,dim(Data)[1])), Data[,-1], stringsAsFactors=F)
outfile = paste0(comp,".datamatrix.txt")
cat("Protein\tLength\t0\t0\t0\t1\t1\t1\n", file=outfile, sep="", append=F)
for (i in 1:dim(D)[1]) {
	cat(unlist(D[i,]), file=outfile, sep="\t", append=T)
	cat("\n", file=outfile, append=T)
}

# // WIN32 VERSION //
# run these commands in the Cygwin terminal window
# (working dir: MISC/R/omics/2015-04-15)
# ../qprot-win32/qprot_1.3.0/qspec-param.exe prep-wt.datamatrix.txt 5000 20000 0
# ../qprot-win32/qprot_1.3.0/getfdr.exe prep-wt.datamatrix.txt_qspec

# clean up the output (remove the "dummy" length column and sort in FDR order)
R = read.table(paste0(comp,".datamatrix.txt_qspec_fdr"), header=T, sep="\t")
R = R[,-2]
names(R)[2:7] = names(Data)[2:7]
R = R[order(R$fdr),]
write.csv(R, file=paste0(comp,".qspec_results.csv"), row.names=FALSE)

