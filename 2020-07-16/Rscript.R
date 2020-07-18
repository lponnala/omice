
# ----------------------------------------
# // Heirarchical Clustering //
# see script.py for preparing the data
# [the following analysis is based on 2014-11-24/analyz.R]

library(readr)
library(gplots)

for (set in c('set1','set2','set3')) {
    for (typ in c('byProtein','byTissue')) {
        cat(paste0("-- ",set,",",typ," --"))
        data_file = paste0(set,'_data_',typ,'.csv')
        dendro_file = paste0(set,'_clusters_',typ,'.png')
        heatmap_file = paste0(set,'_heatmap_',typ,'.png')

        DATA = readr::read_csv(data_file)
        dim(DATA)
        colnames(DATA)
        D = DATA[,-1]
        # sapply(D,class)
        unique(sapply(D,class))

        # ~~ Dendrogram ~~
        save_fig = TRUE
        D_dd = as.dist((1-cor(t(D)))/2)
        D_hc = hclust(D_dd, method="average")
        if (save_fig) {
            png(filename = dendro_file, width=960, height=480, units="px")
        }
        plot(D_hc, labels=unlist(DATA[,1],use.names=FALSE), hang=-1, frame.plot=FALSE, main="Heirarchical Clusters", sub="", xlab="", ylab="correlation-based distance")
        # num_clust = c(set1 = 8, set2 = 8, set3 = 5)[set]
        # rect.hclust(D_hc, k=num_clust, border = 1 + 1:num_clust)
        if (save_fig) {
            dev.off()
        }

        # ~~ Heatmap ~~
        Dm = as.matrix(D)
        png(filename = heatmap_file, width=960, height=480, units="px")
        heatmap.2(Dm, Rowv = as.dendrogram(D_hc), dendrogram = "row", col=redgreen(75), scale="row", key=TRUE, density.info="none", trace="none", cexCol=0.8, labRow=NA)
        dev.off()
    }
}
