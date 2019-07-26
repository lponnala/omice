
# # ==== PCA ====

# library(xlsx)

# FILE = 'PCA-Lalit.xlsx'
# R = read.xlsx(FILE, sheetIndex=1, colIndex=c(1,7:10,12:15), rowIndex=1:227, header=TRUE, as.data.frame=TRUE)
# names(R) = c('BestModel', 'WT-FL-rep1','WT-FL-rep2','met1-FL-rep1','met1-FL-rep2', 'WT-NL-rep1','WT-NL-rep2','met1-NL-rep1','met1-NL-rep2')
# D = R[,-1]

# y = princomp(D)
# print(summary(y))
# print(y$loadings, cutoff=0.0001)
# jpeg('pcaPlot.jpg')
# plot(y$loadings[1:2,1:2],type='p',pch=15,main="PCA using Covariance matrix",xlab='PC1',ylab='PC2',xlim=c(min(y$loadings[,1])-0.1,max(y$loadings[,1])+0.1),ylim=c(min(y$loadings[,2])-0.1,max(y$loadings[,2])+0.1))
# points(y$loadings[3:4,1:2],type='p',pch=17)
# points(y$loadings[5:6,1:2],type='p',pch=2)
# points(y$loadings[7:8,1:2],type='p',pch=7)
# legend('topleft',inset=0.05,legend=c('Wt-FL','met1-FL','Wt-NL','met1-NL'),pch=c(15,17,2,7),col='black')
# dev.off()

# fo = file('pcaCoords.txt','w')
# for (i in 1:dim(y$loadings)[1]) {
	# cat(rownames(y$loadings)[i],'\t',y$loadings[i,1],'\t',y$loadings[i,2],'\n',sep='',file=fo)
# }
# close(fo)

# # y = princomp(D, cor=TRUE)
# # print(summary(y))
# # print(y$loadings, cutoff=0.0001)
# # jpeg('pcaPlot_usingCor.jpg')
# # plot(y$loadings[1:2,1:2],type='p',pch=15,main="PCA using Correlation matrix",xlab='PC1',ylab='PC2',xlim=c(min(y$loadings[,1])-0.1,max(y$loadings[,1])+0.1),ylim=c(min(y$loadings[,2])-0.1,max(y$loadings[,2])+0.1))
# # points(y$loadings[3:4,1:2],type='p',pch=17)
# # points(y$loadings[5:6,1:2],type='p',pch=2)
# # points(y$loadings[7:8,1:2],type='p',pch=7)
# # legend('topleft',inset=0.05,legend=c('Wt-FL','met1-FL','Wt-NL','met1-NL'),pch=c(15,17,2,7),col='black')
# # dev.off()

# # fo = file('pcaCoords_usingCor.txt','w')
# # for (i in 1:dim(y$loadings)[1]) {
	# # cat(rownames(y$loadings)[i],'\t',y$loadings[i,1],'\t',y$loadings[i,2],'\n',sep='',file=fo)
# # }
# # close(fo)


# # ==== Differential expression ====

# source("glee-r-funcs.r")

# # -- user-defined --
# name = "FL" # "NL"
# row_ind = 1:226 # 1:212

# XLFILE = paste(name,".spc-stats-lalit.xls",sep="")

# # number of replicates
# nA = 2
# nB = 2

# # specify rowIndex and/or colIndex if errors occur
# Data = read.xlsx(XLFILE, sheetIndex=1, rowIndex=row_ind, colIndex=NULL, as.data.frame=TRUE, header=TRUE)

# # options
# fit_type = "cubic"
# num_iter = 10000
# num_digits = 4
# outfile = paste(name,".diff-exp.xlsx",sep="")
# # ------------------

# if (!data_ok(Data,nA,nB)) {
	# msg = paste("check input data! spreadsheet must have",
			# "(1) the right number of columns",
			# "(2) positive finite values",
			# "\n",sep="\n")
	# stop(msg)
# }

# Prot = as.character(Data[,1])
# A = as.matrix(Data[,1+(1:nA)])
# B = as.matrix(Data[,1+nA+(1:nB)])

# m = fit_model(A, B, fit_type)
# model_fit_plots(m, outfile="fitplots.png")
# stn_pval = calc_stn_pval(A, B, m, num_iter)
# stn_pval_plots(stn_pval, outfile="stn-pval.png")
# diff_exp_table(stn_pval, Prot, num_digits, outfile)


# ==== Hierarchical clustering ====

library(xlsx)

XLSXFILE = "cluster-lalit-NadjSPC.xlsx"
DATAFILE = "cluster-lalit-NadjSPC.RData"

# colNames = unlist(read.xlsx(XLSXFILE, sheetIndex=1, rowIndex=1, colIndex=1:9, as.data.frame=FALSE, header=FALSE))
# DATA = read.xlsx(XLSXFILE, sheetIndex=1, rowIndex=2:150, colIndex=1:9, as.data.frame=TRUE, header=FALSE)
# names(DATA) = colNames
# save(DATA, file=DATAFILE)

load(DATAFILE)

D = DATA[,6:9]
D_dd = as.dist((1-cor(t(D)))/2)
D_hc = hclust(D_dd, method="average")
# png(filename = "clusters.png", width=960, height=480, units="px")
plot(D_hc, labels=FALSE, hang=-1, frame.plot=FALSE, main="Heirarchical Clusters", sub="", xlab="", ylab="correlation-based distance")
num_clust = 6
rect.hclust(D_hc, k=num_clust, border = 1 + 1:num_clust)
# dev.off()

# # to identify the contents of each cluster
# x = identify(D_hc)
# # click on each cluster manually and then save the resulting list
# save(x, file="clusters.RData")
# # print the contents of each cluster
# for (i in 1:length(x)) {
	# cat("\n", "-- cluster ", i, " --", "\n", sep="")
	# cat(as.character(DATA[x[[i]],1]), sep="\n")
	# cat("\n\n")
# }


# to draw a heatmap
library(gplots)
Dm = as.matrix(D)
colnames(Dm) = c("WT-FL", "Met1-FL", "WT-NL", "Met1-NL")

# helpful links for heatmap
# https://www.biostars.org/p/14156/
# http://mannheimiagoesprogramming.blogspot.com/2012/06/drawing-heatmaps-in-r-with-heatmap2.html

heatmap.2(Dm, Rowv = as.dendrogram(D_hc), dendrogram = "row", col=redgreen(75), scale="row", key=FALSE, density.info="none", trace="none", cexCol=0.8, labRow=NA)

# [works] heatmap(Dm, cexCol=0.7, labRow=NA)
# [works] heatmap.2(Dm, Rowv = as.dendrogram(D_hc))
# heatmap.2(Dm, Rowv = as.dendrogram(D_hc), col=redgreen(75), scale="row", key=FALSE, density.info="none", trace="none", cexCol=0.7, labRow=NA)
# heatmap.2(Dm, col=redgreen(75), scale="row", key=T, keysize=0.5, density.info="none", trace="none", cexCol=0.7, labRow=NA)
# heatmap.2(Dm, col=redgreen(75), dendrogram="row",
          # scale="row", key=T, keysize=0.5, density.info="none",
          # trace="none",cexCol=0.5, labRow=NA)

# cross-check the color coding
load("clusters.RData")
# compare "Met1-FL" between clusters 2 and 5:
cat("Average level for Met1-FL in cluster2 (mostly red) = ", mean(x[[2]],2), "\n", sep="")
cat("Average level for Met1-FL in cluster5 (mostly green) = ", mean(x[[5]],2), "\n", sep="")
