library(xlsx)

# --------------------------------------------------
# WT-clpp3 dataset
# --------------------------------------------------
FILE='PCA-set-wt-clpp3-1.xls'
R = read.xlsx(FILE, sheetIndex=1, rowIndex=seq(1,2117), colIndex=seq(1,7), header=TRUE, as.data.frame=TRUE)
D = R[,-1]
y = princomp(D)
print(summary(y))

print(y$loadings)
plot(y$loadings[1:3,1:2],type='p',pch=15,xlab='PC1',ylab='PC2',xlim=c(-0.6,-0.2),ylim=c(-0.5,0.6))
points(y$loadings[4:6,1:2],type='p',pch=17)
legend('topleft',inset=0.05,legend=c('WT','clpp3'),pch=c(15,17),col='black')

fo = file('pca_coords.txt','w')
for (i in 1:dim(y$loadings)[1]) {
	cat(rownames(y$loadings)[i],'\t',y$loadings[i,1],'\t',y$loadings[i,2],'\n',sep='',file=fo)
}
close(fo)
# --------------------------------------------------


# # --------------------------------------------------
# # OTHER 2 SAMPLES
# # --------------------------------------------------

# FILE = 'correlation-wt-clpcsdouble.xls'
# R = read.xlsx(FILE, sheetIndex=1, rowIndex=seq(1,1078), colIndex=seq(1,7), header=TRUE, as.data.frame=TRUE)
# D = R[,-1]
# y = princomp(D)
# print(summary(y))

# # Importance of components:
	# # Comp.1	Comp.2	Comp.3	Comp.4	Comp.5	Comp.6
# # Standard deviation	0.009267	0.001501	0.00123	0.0007872	0.000711	0.0005533
# # Proportion of Variance	0.942923	0.024741	0.01662	0.0068027	0.005550	0.0033605
# # Cumulative Proportion	0.942923 	0.967665	0.98429	0.9910894	0.996640	1.0000000

# print(y$loadings)
# plot(y$loadings[c(1,3,5),1:2],type='p',pch=16,xlab='PC1',ylab='PC2',xlim=c(-0.48,-0.32),ylim=c(-0.7,0.5))
# points(y$loadings[c(2,4,6),1:2],type='p',pch=17)
# legend('topright',inset=0.05,legend=c('WT','clpC1xS'),pch=c(16,17),col='black')

# fo = file('pca2.txt','w')
# for (i in 1:dim(y$loadings)[1]) {
	# cat(rownames(y$loadings)[i],'\t',y$loadings[i,1],'\t',y$loadings[i,2],'\n',sep='',file=fo)
# }
# close(fo)
# # --------------------------------------------------


# # --------------------------------------------------
# # ALL 3 SAMPLES
# # --------------------------------------------------

# FILE='correlation-wt-clps-clpc.xls'
# R = read.xlsx(FILE, sheetIndex=1, rowIndex=seq(1,791), colIndex=seq(1,10), header=TRUE, as.data.frame=TRUE)
# D = R[,-1]
# y = princomp(D)
# print(summary(y))

# # Importance of components:
	# # Comp.1	Comp.2	Comp.3	Comp.4	Comp.5	Comp.6	Comp.7	Comp.8	Comp.9
# # Standard deviation	0.01283	0.001695	0.001661	0.0009141	0.0006781	0.0006043	 0.0005191	0.0003841	0.0003773
# # Proportion of Variance	0.95446	0.016664	0.016000	0.0048466	0.0026676	0.0021182 	0.0015634	0.0008557	0.0008256
# # Cumulative Proportion	 0.95446	0.971123	0.987123	0.9919697	0.9946372	0.9967554	0.9983188	0.9991744	1.0000000

# print(y$loadings)
# plot(y$loadings[1:3,1:2],type='p',pch=15,xlab='PC1',ylab='PC2',xlim=c(-0.4,-0.25),ylim=c(-0.38,0.52))
# points(y$loadings[4:6,1:2],type='p',pch=16)
# points(y$loadings[7:9,1:2],type='p',pch=17)
# legend('topleft',inset=0.05,legend=c('WT','clpS','clpC'),pch=c(15,16,17),col='black')
	
# fo = file('pca.txt','w')
# for (i in 1:dim(y$loadings)[1]) {
	# cat(rownames(y$loadings)[i],'\t',y$loadings[i,1],'\t',y$loadings[i,2],'\n',sep='',file=fo)
# }
# close(fo)
# # --------------------------------------------------
