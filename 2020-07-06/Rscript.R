
# ----------------------------------------
# // Heirarchical Clustering //
# see script.py for preparing the data
# [the following analysis is based on 2014-11-24/analyz.R]

DATA = readr::read_csv("cluster_data.csv")
dim(DATA)
colnames(DATA)

D = DATA[,-1]
sapply(D,class)
unique(sapply(D,class))

# ~~ Dendrogram ~~
save_fig = TRUE

D_dd = as.dist((1-cor(t(D)))/2)
D_hc = hclust(D_dd, method="average")
if (save_fig) {
    png(filename = "clusters.png", width=960, height=480, units="px")
}
plot(D_hc, labels=unlist(DATA[,1],use.names=FALSE), hang=-1, frame.plot=FALSE, main="Heirarchical Clusters", sub="", xlab="", ylab="correlation-based distance")
num_clust = 4
rect.hclust(D_hc, k=num_clust, border = 1 + 1:num_clust)
if (save_fig) {
    dev.off()
}

# ~~ Heatmap ~~
library(gplots)

Dm = as.matrix(D)
# colnames(Dm) = c("WT-FL", "Met1-FL", "WT-NL", "Met1-NL")

# helpful links for heatmap
# https://www.biostars.org/p/14156/
# http://mannheimiagoesprogramming.blogspot.com/2012/06/drawing-heatmaps-in-r-with-heatmap2.html

png(filename = "heatmap.png", width=960, height=480, units="px")
heatmap.2(Dm, Rowv = as.dendrogram(D_hc), dendrogram = "row", col=redgreen(75), scale="row", key=FALSE, density.info="none", trace="none", cexCol=0.8, labRow=NA)
dev.off()


# ----------------------------------------
# // Differential Expression //
# see script.py for preparing the data
# [the following analysis is based on 2018-02-17/compare-across-funcs.R]

y = readr::read_csv("diffexp_anova_data.csv")

# ----------------------------------------
