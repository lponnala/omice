
# --------------------------------------
# run PCA
# --------------------------------------

# -- using adjSPC --
load("A.RData")
y = princomp(A[,-1])
print(summary(y))
print(y$loadings, cutoff=0.0001)
jpeg('adjSPC-pcaPlot.jpg')
plot(y$loadings[1:3,1:2],type='p',pch=15,main="PCA using adjSPC data",xlab='PC1',ylab='PC2',xlim=c(min(y$loadings[,1])-0.1,max(y$loadings[,1])+0.1),ylim=c(min(y$loadings[,2])-0.1,max(y$loadings[,2])+0.1))
points(y$loadings[4:6,1:2],type='p',pch=17)
legend('topleft',inset=0.05,legend=c('wt','p1xp2'),pch=c(15,17),col='black')
dev.off()
write.csv(y$loadings[,1:2], file = "adjSPC-pca.csv")

# -- using NadjSPC --
load("N.RData")
y = princomp(N[,-1])
print(summary(y))
print(y$loadings, cutoff=0.0001)
jpeg('NadjSPC-pcaPlot.jpg')
plot(y$loadings[1:3,1:2],type='p',pch=15,main="PCA using NadjSPC data",xlab='PC1',ylab='PC2',xlim=c(min(y$loadings[,1])-0.1,max(y$loadings[,1])+0.1),ylim=c(min(y$loadings[,2])-0.1,max(y$loadings[,2])+0.1))
points(y$loadings[4:6,1:2],type='p',pch=17)
legend('topleft',inset=0.05,legend=c('wt','p1xp2'),pch=c(15,17),col='black')
dev.off()
write.csv(y$loadings[,1:2], file = "NadjSPC-pca.csv")

# --------------------------------------


# # --------------------------------------
# # run QSPEC
# # --------------------------------------

# # using adjSPC
# load("A.RData")
# Data = A
# comp = "adjSPC"

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

# # // WIN32 VERSION //
# # run these commands in the Cygwin terminal window
# # (working dir: MISC/R/omics/2015-07-24)
# # ../qprot-win32/qprot_1.3.0/qspec-param.exe adjSPC.datamatrix.txt 5000 20000 0
# # ../qprot-win32/qprot_1.3.0/getfdr.exe adjSPC.datamatrix.txt_qspec

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
# load("A.RData")

# # -- using shiny app --
# # run using 5000 iterations at: http://lponnala.shinyapps.io/glee
# write.table(A, file="A.txt", quote=FALSE, sep="\t", eol = "\n", row.names=FALSE)
# # show all differentially-expressed proteins, copy-paste the html output into textfile, then copy-paste text into spreadsheet

# # # -- using gleeR package --
# # library(gleeR)
# # # [failed to run due to memory outage]
# # # options
# # Data = A
# # comp = "adjSPC"
# # nA = 3
# # nB = 3
# # fit_type = "cubic"
# # num_iter = 1000
# # num_digits = 4
# # outfile = ""
# # # run
# # Prot = as.character(Data[,1])
# # A = as.matrix(Data[,1+(1:nA)])
# # B = as.matrix(Data[,1+nA+(1:nB)])
# # if (!data_ok(Data,nA,nB)) {
	# # msg = paste("check input data! spreadsheet must have",
			# # "(1) the right number of columns",
			# # "(2) positive finite values",
			# # "\n",sep="\n")
	# # stop(msg)
# # }
# # m = fit_model(A, B, fit_type)
# # model_fit_plots(m, outfile=paste0(comp,".fitplots.png"))
# # stn_pval = calc_stn_pval(A, B, m, num_iter)
# # stn_pval_plots(stn_pval, outfile=paste0(comp,".stn-pval.png"))
# # T = diff_exp_table(stn_pval, Prot, num_digits, outfile)
# # write.csv(T, file=paste0(comp,".glee-diffexp.csv"), row.names=FALSE)

# # --------------------------------------


# # --------------------------------------
# # run correlations
# # --------------------------------------
# load("A.RData")
# write.csv(cor(A[,-1]), file="corr-adjSPC.csv")
# load("N.RData")
# write.csv(cor(N[,-1]), file="corr-NadjSPC.csv")
# # --------------------------------------


# # --------------------------------------
# # prepare the data
# # --------------------------------------
# library(readxl)
# D = read_excel("forlalit-wt-vsprepdouble.xlsx", sheet = 1, col_names = TRUE)
# D = D[1:913,c(1:2,4:15)]
# A = D[,c(2,6:8,12:14)]
# names(A) = c("Acc","wt-rep1","wt-rep2","wt-rep3","p1xp2-rep1","p1xp2-rep2","p1xp2-rep3")
# save(A,file="A.RData")
# N = D[,c(2,3:5,9:11)]
# names(N) = c("Acc","wt-rep1","wt-rep2","wt-rep3","p1xp2-rep1","p1xp2-rep2","p1xp2-rep3")
# save(N,file="N.RData")
# # --------------------------------------

