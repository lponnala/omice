
# ----------------------------------------
# // Heirarchical Clustering //
# see script.py for preparing the data

# heatmap: check out pheatmap()
# https://www.datanovia.com/en/blog/clustering-using-correlation-as-distance-measures-in-r/


# -- cluster the proteins (show p-values) --
library(readr)
library(pvclust)
set.seed(123)
num_boot = 1000
save_plots = TRUE
for (key in c("PGs","17-ABC1Ks")) {
    for (typ in c("Scaled","Zscore")) {
        cat(paste0("\n","-- ",key,",",typ," --","\n"))

        # ~~ Data ~~
        data_file = paste0(key,'_',typ,'_data.csv')
        DATA = readr::read_csv(data_file)
        D = t(DATA[,-1])
        colnames(D) = unlist(DATA[,1])
        unique(sapply(D,class))

        # ~~ Dendrogram ~~
        # pv = pvclust(D, method.hclust="average", method.dist="correlation", nboot=20) # use the default correlation-distance (which does not divide by 2)
        pv = pvclust(D, method.hclust="average", method.dist = function(x) as.dist((1-cor(x))/2), nboot=num_boot) # divide by 2 to make dist go from 0 to 1

        # ~~ Plot ~~
        if (save_plots) {
            dendro_file = paste0(key,'_',typ,'_clusters_pval.png')
            png(filename = dendro_file, width=960, height=480, units="px")
        }
        plot(pv, print.pv="au", print.num=FALSE, hang = -1, col.pv=c(au="blue"), cex.pv=1, cex = 1, frame.plot=FALSE, main=paste0("Heirarchical Clusters: ",key," (using ",typ," abundance)"), sub="", xlab="", ylab="correlation-based distance")
        pvrect(pv, alpha=0.95, pv="au") # mark red boxes around the significant clusters
        if (save_plots) {
            dev.off()
        }
    }
}


# # -- cluster the proteins (no p-values) --
# library(readr)
# for (key in c("PGs","17-ABC1Ks")) {
#     for (typ in c("Scaled","Zscore")) {
#         cat(paste0("\n","-- ",key,",",typ," --","\n"))
#         data_file = paste0(key,'_',typ,'_data.csv')
#         dendro_file = paste0(key,'_',typ,'_clusters.png')
# 
#         DATA = readr::read_csv(data_file)
#         dim(DATA)
#         colnames(DATA)
#         D = DATA[,-1]
#         # sapply(D,class)
#         unique(sapply(D,class))
#         # print(head(D))
# 
#         # ~~ Dendrogram ~~
#         D_dd = as.dist((1-cor(t(D)))/2) # divide by 2 so that the max dist is 1
#         # D_dd = as.dist(1-cor(t(D))) # leave out the "dividing by 2" part
#         # print(D_dd)
#         D_hc = hclust(D_dd, method="average")
#         png(filename = dendro_file, width=960, height=480, units="px")
#         plot(D_hc, labels=unlist(DATA[,1],use.names=FALSE), hang=-1, frame.plot=FALSE, main=paste0("Heirarchical Clusters: ",key," (using ",typ," abundance)"), sub="", xlab="", ylab="correlation-based distance")
#         # num_clust = 4
#         # rect.hclust(D_hc, k=num_clust, border = 1 + 1:num_clust)
#         dev.off()
#     }
# }


# -- cluster the tissues (no p-values) --
library(readr)
for (key in c("PGs","17-ABC1Ks")) {
    for (typ in c("Scaled","Zscore")) {
        cat(paste0("\n","-- ",key,",",typ," --","\n"))
        data_file = paste0(key,'_',typ,'_data.csv')
        dendro_file = paste0(key,'_',typ,'_tissue_clusters.png')

        DATA = readr::read_csv(data_file)
        dim(DATA)
        colnames(DATA)
        D = t(DATA[,-1])
        colnames(D) = unlist(DATA[,1])

        # ~~ Dendrogram ~~
        D_dd = as.dist((1-cor(t(D)))/2) # divide by 2 so that the max dist is 1
        # D_dd = as.dist(1-cor(t(D))) # leave out the "dividing by 2" part
        # print(D_dd)
        D_hc = hclust(D_dd, method="average")
        png(filename = dendro_file, width=960, height=480, units="px")
        plot(D_hc, labels=colnames(DATA)[-1], hang=-1, frame.plot=FALSE, main=paste0("Heirarchical Clusters: ",key," (using ",typ," abundance)"), sub="", xlab="", ylab="correlation-based distance")
        # num_clust = 4
        # rect.hclust(D_hc, k=num_clust, border = 1 + 1:num_clust)
        dev.off()
    }
}


# # -- heatmap --
# library(gplots)
# key = 'PGs'
# typ = 'Scaled'
# data_file = paste0(key,'_',typ,'_data.csv')
# DATA = readr::read_csv(data_file)
# D = DATA[,-1]
# hclust_prot = hclust(as.dist((1-cor(t(D)))/2), method="average")
# hclust_tiss = hclust(as.dist((1-cor(D))/2), method="average")
# heatmap.2(as.matrix(D), scale = "row", col = bluered(100), 
#           trace = "none", density.info = "none",
#           Colv = as.dendrogram(hclust_tiss),
#           Rowv = as.dendrogram(hclust_prot))

