library(xlsx)

FILE = 'raw-data.xls' # 'forlalit-ClpT1T1-wtnewtair10.xlsx'
R = read.xlsx(FILE, sheetIndex=1, colIndex=c(1,2:7), rowIndex=1:1882, header=TRUE, as.data.frame=TRUE)
names(R) = c('accession','WT_1','WT_2','WT_3','ClpT_1','ClpT_2','ClpT_3')
D = R[1:1881,-1]


# # ------------------------------------------
# # // GLEE //
# # ------------------------------------------
# # first, run GLEE and then the following steps
# library(fdrtool)
# P = read.xlsx('./glee_output/10k/GLEE_DEG.xlsx', sheetIndex=1, header=TRUE, colIndex=13)
# y = fdrtool(P[,1], statistic="pvalue", plot=FALSE)
# write.xlsx(cbind(y$pval,y$qval), file='out.xlsx')
# # ------------------------------------------


# ------------------------------------------
# // PCA //
# ------------------------------------------
y = princomp(D)
print(summary(y))
print(y$loadings, cutoff=0.001)
jpeg('pcaPlot.jpg')
plot(y$loadings[1:3,1:2],type='p',pch=15,xlab='PC1',ylab='PC2',xlim=c(min(y$loadings[,1])-0.1,max(y$loadings[,1])+0.1),ylim=c(min(y$loadings[,2])-0.1,max(y$loadings[,2])+0.1))
points(y$loadings[4:6,1:2],type='p',pch=17)
legend('topleft',inset=0.05,legend=c('WT','ClpT'),pch=c(15,17),col='black')
dev.off()

fo = file('pcaCoords.txt','w')
for (i in 1:dim(y$loadings)[1]) {
	cat(rownames(y$loadings)[i],'\t',y$loadings[i,1],'\t',y$loadings[i,2],'\n',sep='',file=fo)
}
close(fo)

y = princomp(D, cor=TRUE)
print(summary(y))
print(y$loadings)
jpeg('pcaPlot_usingCor.jpg')
plot(y$loadings[1:3,1:2],type='p',pch=15,xlab='PC1',ylab='PC2',xlim=c(min(y$loadings[,1])-0.1,max(y$loadings[,1])+0.1),ylim=c(min(y$loadings[,2])-0.1,max(y$loadings[,2])+0.1))
points(y$loadings[4:6,1:2],type='p',pch=17)
legend('topleft',inset=0.05,legend=c('WT','ClpT'),pch=c(15,17),col='black')
dev.off()

fo = file('pcaCoords_usingCor.txt','w')
for (i in 1:dim(y$loadings)[1]) {
	cat(rownames(y$loadings)[i],'\t',y$loadings[i,1],'\t',y$loadings[i,2],'\n',sep='',file=fo)
}
close(fo)
# ------------------------------------------


# ------------------------------------------
# // CORRELATIONS //
# ------------------------------------------

fo = file('corr.txt','w')

cat('\n\n\n---- BETWEEN REPLICATES ----\n',sep='', file=fo);

# ---- WT ----
X = D[,1:3]
cat('\n\nWT\tNUM_PROTEINS\tCORRELATION\tPVALUE\n', file=fo)
i=1; j=2;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat('Rep',i,',Rep',j,'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)
# print('Rep%d,Rep%d\t%d\t%f\t%f\n',i,j,length(I),y[[4]],y[[3]])
i=2; j=3;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat('Rep',i,',Rep',j,'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)
i=1; j=3;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat('Rep',i,',Rep',j,'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)

# ---- ClpT ----
X = D[,4:6]
cat('\n\nClpT\tNUM_PROTEINS\tCORRELATION\tPVALUE\n', file=fo)
i=1; j=2;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat('Rep',i,',Rep',j,'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)
i=2; j=3;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat('Rep',i,',Rep',j,'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)
i=1; j=3;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat('Rep',i,',Rep',j,'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)

cat('\n\n\n---- BETWEEN GENOTYPES ----\n',sep='', file=fo);
G=c('WT','ClpT')

# ---- Rep1 ----
X=D[,c(1,4)]
cat('\n\nRep1\tNUM_PROTEINS\tCORRELATION\tPVALUE\n', file=fo)
i=1; j=2;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat(G[i],',',G[j],'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)

# ---- Rep2 ----
X=D[,c(2,5)]
cat('\n\nRep2\tNUM_PROTEINS\tCORRELATION\tPVALUE\n', file=fo)
i=1; j=2;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat(G[i],',',G[j],'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)

# ---- Rep3 ----
X=D[,c(3,6)]
cat('\n\nRep3\tNUM_PROTEINS\tCORRELATION\tPVALUE\n', file=fo)
i=1; j=2;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat(G[i],',',G[j],'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)

close(fo)
# ------------------------------------------
