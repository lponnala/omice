
# --------------------------------------------------
# OTHER 2 SAMPLES
# --------------------------------------------------

FILE = 'correlation-wt-clpcsdouble.xls'
R = read.xlsx(FILE, sheetIndex=1, rowIndex=seq(1,1078), colIndex=seq(1,7), header=TRUE, as.data.frame=TRUE)
D = R[,-1]

fo = file('corr2.txt','w')

cat('\n\n\n---- BETWEEN REPLICATES ----\n',sep='', file=fo);

# ---- WT ----
X = D[,c(1,3,5)]
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

# ---- clpC1xS ----
X = D[,c(2,4,6)]
cat('\n\nclpC1xS\tNUM_PROTEINS\tCORRELATION\tPVALUE\n', file=fo)
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
G=c('WT','clpC1xS')

# ---- Rep1 ----
X=D[,c(1,2)]
cat('\n\nRep1\tNUM_PROTEINS\tCORRELATION\tPVALUE\n', file=fo)
i=1; j=2;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat(G[i],',',G[j],'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)

# ---- Rep2 ----
X=D[,c(3,4)]
cat('\n\nRep2\tNUM_PROTEINS\tCORRELATION\tPVALUE\n', file=fo)
i=1; j=2;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat(G[i],',',G[j],'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)

# ---- Rep3 ----
X=D[,c(5,6)]
cat('\n\nRep3\tNUM_PROTEINS\tCORRELATION\tPVALUE\n', file=fo)
i=1; j=2;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat(G[i],',',G[j],'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)

close(fo)
# --------------------------------------------------


# --------------------------------------------------
# ALL 3 SAMPLES
# --------------------------------------------------

FILE='correlation-wt-clps-clpc.xls'
R = read.xlsx(FILE, sheetIndex=1, rowIndex=seq(1,791), colIndex=seq(1,10), header=TRUE, as.data.frame=TRUE)
D = R[,-1]

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

# ---- clpS ----
X = D[,4:6]
cat('\n\nclpS\tNUM_PROTEINS\tCORRELATION\tPVALUE\n', file=fo)
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

# ---- clpC ----
X = D[,7:9]
cat('\n\nclpC\tNUM_PROTEINS\tCORRELATION\tPVALUE\n', file=fo)
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
G=c('WT','clpS','clpC')

# ---- Rep1 ----
X=D[,c(1,4,7)]
cat('\n\nRep1\tNUM_PROTEINS\tCORRELATION\tPVALUE\n', file=fo)
i=1; j=2;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat(G[i],',',G[j],'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)
i=2; j=3;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat(G[i],',',G[j],'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)
i=1; j=3;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat(G[i],',',G[j],'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)

# ---- Rep2 ----
X=D[,c(2,5,8)]
cat('\n\nRep2\tNUM_PROTEINS\tCORRELATION\tPVALUE\n', file=fo)
i=1; j=2;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat(G[i],',',G[j],'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)
i=2; j=3;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat(G[i],',',G[j],'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)
i=1; j=3;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat(G[i],',',G[j],'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)

# ---- Rep3 ----
X=D[,c(3,6,9)]
cat('\n\nRep3\tNUM_PROTEINS\tCORRELATION\tPVALUE\n', file=fo)
i=1; j=2;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat(G[i],',',G[j],'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)
i=2; j=3;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat(G[i],',',G[j],'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)
i=1; j=3;
I = which(X[,i]+X[,j] > 0)
y = cor.test(X[I,i],X[I,j])
cat(G[i],',',G[j],'\t',length(I),'\t',sprintf('%1.6f',y[[4]]),'\t',sprintf('%1.6f',y[[3]]),'\n', sep='',file=fo)

close(fo)
# --------------------------------------------------
