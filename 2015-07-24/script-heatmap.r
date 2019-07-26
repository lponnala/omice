
# --------------------------------------
# draw heatmap with/without dendrogram
# --------------------------------------

# library(readxl)
# D = read_excel("heatmap.xlsx", sheet = 1, col_names = TRUE)
# D = D[1:50,]
# save(D, file="HeatmapData.RData")

library(gplots)
load("HeatmapData.RData")

A = as.matrix(D[,2:5])
colnames(A) = c("wt","prep","oop","triple")
rownames(A) = as.character(D[,1])
A_dd = as.dist((1-cor(t(A)))/2)
A_hc = hclust(A_dd, method="average")
heatmap.2(A, Rowv = as.dendrogram(A_hc), dendrogram = "row", col=redgreen(75), scale="row", key=FALSE, margins=c(5,8), density.info="none", trace="none")

# divide into 9 parts (should match with the marking sent by Klaas)
plot(A_hc, labels=FALSE, hang=-1, frame.plot=FALSE, main="Heirarchical Clusters", sub="", xlab="", ylab="correlation-based distance")
num_clust = 9
rect.hclust(A_hc, k=num_clust, border = 1 + 1:num_clust)

# to identify the contents of each cluster
x = identify(A_hc)
# click on each cluster manually and then save the resulting list
save(x, file="HeatmapClusters.RData")
# print the contents of each cluster
for (i in 1:length(x)) {
	cat("\n", "-- cluster ", i, " --", "\n", sep="")
	# cat(as.character(A[x[[i]],1]), sep="\n")
	cat(rownames(A)[x[[i]]], sep="\n")
	cat("\n\n")
}

# S = as.matrix(D[,6:9])
# colnames(S) = c("wt","prep","oop","triple")
# rownames(S) = as.character(D[,1])
# heatmap.2(S, Rowv=NULL, Colv=NULL, col=bluered(255), dendrogram="none", scale="none", trace="none", key=TRUE, margins=c(5,8), colsep=1:ncol(S), sepcolor="black", sepwidth=0.005, main="Stdev NadjSPC")

# C = as.matrix(D[,10:13])
# colnames(C) = c("wt","prep","oop","triple")
# rownames(C) = as.character(D[,1])
# heatmap.2(C, Rowv=NULL, Colv=NULL, col=bluered(255), dendrogram="none", scale="none", trace="none", key=TRUE, margins=c(5,8), colsep=1:ncol(C), sepcolor="black", sepwidth=0.005, main="CV NadjSPC")

# --------------------------------------

