
# ----------------------------------------
# // Heirarchical Clustering //
# see script.py for preparing the data

library(readr)
library(gplots)

for (key in c("PGs","17-ABC1Ks")) {
    for (typ in c("Scaled","Zscore")) {
        cat(paste0("\n","-- ",key,",",typ," --","\n"))
        data_file = paste0(key,'_',typ,'_data.csv')
        dendro_file = paste0(key,'_',typ,'_clusters.png')

        DATA = readr::read_csv(data_file)
        dim(DATA)
        colnames(DATA)
        D = DATA[,-1]
        # sapply(D,class)
        unique(sapply(D,class))
        print(head(D))

        # ~~ Dendrogram ~~
        D_dd = as.dist((1-cor(t(D)))/2)
        D_hc = hclust(D_dd, method="average")
        png(filename = dendro_file, width=960, height=480, units="px")
        plot(D_hc, labels=unlist(DATA[,1],use.names=FALSE), hang=-1, frame.plot=FALSE, main=paste0("Heirarchical Clusters: ",key," ",typ), sub="", xlab="", ylab="correlation-based distance")
        num_clust = 4
        rect.hclust(D_hc, k=num_clust, border = 1 + 1:num_clust)
        dev.off()
    }
}

