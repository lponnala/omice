
# ----------------------------------------
# // Heirarchical Clustering //
# see script.py for preparing the data

# heatmap: check out pheatmap()
# https://www.datanovia.com/en/blog/clustering-using-correlation-as-distance-measures-in-r/


# -- with p-values --
library(readr)
library(pvclust)
set.seed(123)
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
        # pv = pvclust(D, method.hclust="average", method.dist="correlation", nboot=20)
        pv = pvclust(D, method.hclust="average", method.dist = function(x) as.dist((1-cor(x))/2) , nboot=20)

        # ~~ Plot ~~
        if (save_plots) {
            dendro_file = paste0(key,'_',typ,'_clusters_pval.png')
            png(filename = dendro_file, width=960, height=480, units="px")
        }
        plot(pv, hang = -1, cex = 1, frame.plot=FALSE, main=paste0("Heirarchical Clusters: ",key," (using ",typ," abundance)"), sub="", xlab="", ylab="correlation-based distance")
        pvrect(pv)
        if (save_plots) {
            dev.off()
        }
    }
}


# # -- initial version (no p-values) --
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

