
# ----------------------------------------------
# QSPEC
# ----------------------------------------------

library(xlsx)
FILE = "forlalitUVR-PCA-Stats.xlsx"
R = read.xlsx(FILE, sheetIndex=1, colIndex=c(1,18:26), rowIndex=2:974, header=FALSE, as.data.frame=TRUE)
names(R) = c('accession','wt_1','wt_2','wt_3','UVR_1','UVR_2','UVR_3','ClpS_1','ClpS_2','ClpS_3')

if (!all(is.finite( as.matrix(R[,-1]) ))) { stop("not all NadjSPC values are finite!") }

# comp = "wt-UVR"; D = R[,c(1,2:4,5:7)]
# comp = "wt-ClpS"; D = R[,c(1,2:4,8:10)]
comp = "UVR-ClpS"; D = R[,c(1,5:7,8:10)]

for (j in 1:dim(D)[2]) {
	D[,j] = as.character(D[,j])
}
# replace spaces in protein names
D[,1] = stringr::str_replace_all(D[,1],pattern=" ",replacement="")
# insert a column of ones as a dummy for length
Data = cbind(D[,1], as.character(rep(1,dim(D)[1])), D[,-1], stringsAsFactors=F)
outfile = paste0(comp,".datamatrix.txt")
cat("Protein\tLength\t0\t0\t0\t1\t1\t1\n", file=outfile, sep="", append=F)
for (i in 1:dim(Data)[1]) {
	cat(unlist(Data[i,]), file=outfile, sep="\t", append=T)
	cat("\n", file=outfile, append=T)
}

# // WIN32 VERSION //
# run these commands in the Cygwin terminal window
# (working dir: MISC/omics/2014-10-14)
# ../qprot-win32/qprot_1.3.0/qspec-param.exe wt-UVR.datamatrix.txt 5000 20000 0
# ../qprot-win32/qprot_1.3.0/getfdr.exe wt-UVR.datamatrix.txt_qspec

# clean up the output (remove the "dummy" length column and sort in FDR order)
T = read.table(paste0(comp,".datamatrix.txt_qspec_fdr"), header=T, sep="\t")
T = T[,-2]
names(T)[2:7] = names(D)[2:7]
T = T[order(T$fdr),]
write.csv(T, file=paste0(comp,".qspec_results.csv"), row.names=FALSE)

# ----------------------------------------------


# # ----------------------------------------------
# # GLEE
# # ----------------------------------------------

# library(xlsx)
# FILE = "forlalitUVR-PCA-Stats.xlsx"
# R = read.xlsx(FILE, sheetIndex=1, colIndex=c(1,18:26), rowIndex=2:974, header=FALSE, as.data.frame=TRUE)
# names(R) = c('accession','wt_1','wt_2','wt_3','UVR_1','UVR_2','UVR_3','ClpS_1','ClpS_2','ClpS_3')

# if (!all(is.finite( as.matrix(R[,-1]) ))) { stop("not all NadjSPC values are finite!") }

# D = R[,c(1,2:4,5:7)]
# write.xlsx(D, file="wt-UVR.xlsx", col.names=TRUE, row.names=FALSE)

# D = R[,c(1,2:4,8:10)]
# write.xlsx(D, file="wt-ClpS.xlsx", col.names=TRUE, row.names=FALSE)

# D = R[,c(1,5:7,8:10)]
# write.xlsx(D, file="UVR-ClpS.xlsx", col.names=TRUE, row.names=FALSE)

# # ----------------------------------------------


# # ----------------------------------------------
# # PCA
# # ----------------------------------------------

# # library(xlsx)
# # FILE = "forlalitUVR-PCA-Stats.xlsx"
# # D = read.xlsx(FILE, sheetIndex=1, colIndex=c(1,18:26), rowIndex=2:974, header=FALSE, as.data.frame=TRUE)
# # names(D) = c('accession','wt_1','wt_2','wt_3','UVR_1','UVR_2','UVR_3','ClpS_1','ClpS_2','ClpS_3')
# # save(D, file="dataForPCA.RData")

# load("dataForPCA.RData")

# # calculate principal components
# # # -- using the correlation matrix --
# # y = princomp(D[,-1], cor=TRUE)
# # outname = 'pca-usingCorr'
# # -- using the covariance matrix --
# y = princomp(D[,-1], cor=FALSE)
# outname = 'pca'

# # plot symbols
# pchVals = c(1,8,9)

# # plot PC1 and PC2
# print(summary(y))
# print(y$loadings)
# jpeg(paste(outname,".jpg",sep=""))
# plot(y$loadings[1:3,1:2],type='p',pch=pchVals[1],xlab='PC1',ylab='PC2',xlim=c(min(y$loadings[,1])-0.1,max(y$loadings[,1])+0.1),ylim=c(min(y$loadings[,2])-0.1,max(y$loadings[,2])+0.1))
# points(y$loadings[4:6,1:2],type='p',pch=pchVals[2])
# points(y$loadings[7:9,1:2],type='p',pch=pchVals[3])
# legend('topleft',inset=0.05,legend=c('wt','UVR','ClpS'),pch=pchVals,col='black')
# dev.off()

# # print out PC1 and PC2
# fo = file(paste(outname,".txt",sep=""),'w')
# for (i in 1:dim(y$loadings)[1]) {
	# cat(rownames(y$loadings)[i],'\t',y$loadings[i,1],'\t',y$loadings[i,2],'\n',sep='',file=fo)
# }
# close(fo)

# # ----------------------------------------------


# # ----------------------------------------------
# # Sum the groups and write to Excel spreadsheet
# # ----------------------------------------------

# library(xlsx)
# ExcelFILE = "forlalitUVR.xlsx"
# COLS = unlist(read.xlsx(ExcelFILE, sheetIndex=1, rowIndex=1, colIndex=1:20, as.data.frame=FALSE, header=FALSE))
# DATA = read.xlsx(ExcelFILE, sheetIndex=1, rowIndex=2:126, colIndex=1:20, as.data.frame=TRUE, header=FALSE)
# D = DATA[,c(2,6:8,11:13,16:18)]
# 
# if (!all(sapply(D, FUN = function(x) { all(is.na(x[!is.finite(x)])) } ))) {
# 	stop("there are non-blank cells in the excel spreadsheet that contain non-finite values")
# }
# 
# library(dplyr)
# DF = tbl_df(D)
# G = group_by(D, X2)
# T = summarise(G, count=n(), 
# 							totX6=sum(X6,na.rm=TRUE),
# 							totX7=sum(X7,na.rm=TRUE),
# 							totX8=sum(X8,na.rm=TRUE),
# 							
# 							totX11=sum(X11,na.rm=TRUE),
# 							totX12=sum(X12,na.rm=TRUE),
# 							totX13=sum(X13,na.rm=TRUE),
# 							
# 							totX16=sum(X16,na.rm=TRUE),
# 							totX17=sum(X17,na.rm=TRUE),
# 							totX18=sum(X18,na.rm=TRUE)
# 							)
# names(T) = c("Group","# Proteins",COLS[c(6:8,11:13,16:18)])
# write.xlsx(T, file="summed.xlsx", col.names=TRUE, row.names=FALSE)

# # ----------------------------------------------
