
load("N.RData")

# -- using all data --
y = princomp(N)
print(summary(y))
print(y$loadings, cutoff=0.0001)
jpeg('pcaPlot.jpg')
plot(y$loadings[c(1,5,9),1:2],type='p',pch=15,main="PCA using Covariance matrix",xlab='PC1',ylab='PC2',xlim=c(min(y$loadings[,1])-0.1,max(y$loadings[,1])+0.1),ylim=c(min(y$loadings[,2])-0.1,max(y$loadings[,2])+0.1))
points(y$loadings[c(2,6,10),1:2],type='p',pch=17)
points(y$loadings[c(3,7,11),1:2],type='p',pch=2)
points(y$loadings[c(4,8,12),1:2],type='p',pch=7)
legend('topleft',inset=0.05,legend=c('wt','prep','oop','triple'),pch=c(15,17,2,7),col='black')
dev.off()
write.csv(y$loadings[,1:2], file = "pca.csv")

# -- using totalAdjSPC >= 60 --
load("totalA.RData")
N = N[which(totalA >= 60),]
y = princomp(N)
print(summary(y))
print(y$loadings, cutoff=0.0001)
jpeg('pcaPlot-ge60.jpg')
plot(y$loadings[c(1,5,9),1:2],type='p',pch=15,main="PCA using Covariance matrix",xlab='PC1',ylab='PC2',xlim=c(min(y$loadings[,1])-0.1,max(y$loadings[,1])+0.1),ylim=c(min(y$loadings[,2])-0.1,max(y$loadings[,2])+0.1))
points(y$loadings[c(2,6,10),1:2],type='p',pch=17)
points(y$loadings[c(3,7,11),1:2],type='p',pch=2)
points(y$loadings[c(4,8,12),1:2],type='p',pch=7)
legend('topleft',inset=0.05,legend=c('wt','prep','oop','triple'),pch=c(15,17,2,7),col='black')
dev.off()
write.csv(y$loadings[,1:2], file = "pca-ge60.csv")
