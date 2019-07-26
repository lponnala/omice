
# --------------------------------------
# run PCA
# --------------------------------------

# -- using adjSPC --
load("A_anton.RData")
y = princomp(A[,-1])
print(summary(y))
print(y$loadings, cutoff=0.0001)
jpeg('Anton.adjSPC-pcaPlot.jpg')
plot(y$loadings[1:3,1:2],type='p',pch=15,main="PCA using adjSPC data",xlab='PC1',ylab='PC2',xlim=c(min(y$loadings[,1])-0.1,max(y$loadings[,1])+0.1),ylim=c(min(y$loadings[,2])-0.1,max(y$loadings[,2])+0.1))
points(y$loadings[4:6,1:2],type='p',pch=17)
legend('topleft',inset=0.05,legend=c('wt','mut'),pch=c(15,17),col='black')
dev.off()
write.csv(y$loadings[,1:2], file = "Anton.adjSPC-pca.csv")

# -- using NadjSPC --
load("N_anton.RData")
y = princomp(N[,-1])
print(summary(y))
print(y$loadings, cutoff=0.0001)
jpeg('Anton.NadjSPC-pcaPlot.jpg')
plot(y$loadings[1:3,1:2],type='p',pch=15,main="PCA using NadjSPC data",xlab='PC1',ylab='PC2',xlim=c(min(y$loadings[,1])-0.1,max(y$loadings[,1])+0.1),ylim=c(min(y$loadings[,2])-0.1,max(y$loadings[,2])+0.1))
points(y$loadings[4:6,1:2],type='p',pch=17)
legend('topleft',inset=0.05,legend=c('wt','mut'),pch=c(15,17),col='black')
dev.off()
write.csv(y$loadings[,1:2], file = "Anton.NadjSPC-pca.csv")

# --------------------------------------


# # --------------------------------------
# # run QSPEC
# # --------------------------------------

# # using adjSPC
# load("A_anton.RData")
# Data = A
# comp = "Anton.adjSPC"

# for (j in 1:dim(Data)[2]) {
	# Data[,j] = as.character(Data[,j])
# }
# # insert a column of ones as a dummy for length
# D = cbind(Data[,1], as.character(rep(1,dim(Data)[1])), Data[,-1], stringsAsFactors=FALSE)
# outfile = paste0(comp,".datamatrix.txt")
# cat("Protein\tLength\t0\t0\t0\t1\t1\t1\n", file=outfile, sep="", append=F)
# for (i in 1:dim(D)[1]) {
	# cat(unlist(D[i,]), file=outfile, sep="\t", append=T)
	# cat("\n", file=outfile, append=T)
# }

# # // WIN64 VERSION //
# # run these commands in the Cygwin terminal window
# # (working dir: MISC/R/omics/2015-07-24)
# # ../qprot-win64/qprot_1.3.0/qspec-param.exe Anton.adjSPC.datamatrix.txt 5000 20000 0
# # ../qprot-win64/qprot_1.3.0/getfdr.exe Anton.adjSPC.datamatrix.txt_qspec

# # clean up the output (remove the "dummy" length column and sort in FDR order)
# R = read.table(paste0(comp,".datamatrix.txt_qspec_fdr"), header=T, sep="\t")
# R = R[,-2]
# names(R)[2:7] = names(Data)[2:7]
# R = R[order(R$fdr),]
# write.csv(R, file=paste0(comp,".qspec_results.csv"), row.names=FALSE)

# # --------------------------------------


# # --------------------------------------
# # run GLEE
# # --------------------------------------
# # using adjSPC
# load("A_anton.RData")

# # -- using shiny app --
# # run using 3000 iterations at: http://lponnala.shinyapps.io/glee
# write.table(A, file="A_anton.txt", quote=FALSE, sep="\t", eol = "\n", row.names=FALSE)
# # show all differentially-expressed proteins, copy-paste the html output into textfile, then copy-paste text into spreadsheet
# # --------------------------------------


# # --------------------------------------
# # run correlations
# # --------------------------------------
# load("A_anton.RData")
# write.csv(cor(A[,-1]), file="Anton.corr-adjSPC.csv")
# load("N_anton.RData")
# write.csv(cor(N[,-1]), file="Anton.corr-NadjSPC.csv")
# # --------------------------------------


# # --------------------------------------
# # prepare the data
# # --------------------------------------
# library(readxl)
# D = read_excel("forLalit-Antonset-wt-mutant.xlsx")
# A = D[,c(2,12:17)]
# A[which(is.na(A), arr.ind = TRUE)] = 0
# names(A) = c("Acc","wt-1","wt-2","wt-3","mut-1","mut-2","mut-3")
# save(A, file="A_anton.RData")
# N = D[,c(2,5:10)]
# names(N) = c("Acc","wt-1","wt-2","wt-3","mut-1","mut-2","mut-3")
# save(N, file="N_anton.RData")
# # --------------------------------------
