
# ----------------------------------------
# // Heirarchical Clustering //
# see 2014-11-24/analyz.R

DATA = readr::read_csv("cluster_data.csv")
dim(DATA)
colnames(DATA)

D = DATA[,-1]
sapply(D,class)
unique(sapply(D,class))

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

# # to identify the contents of each cluster
# x = identify(D_hc)
# # click on each cluster manually and then save the resulting list
# save(x, file="clusters.RData")
# # print the contents of each cluster
# for (i in 1:length(x)) {
# 	cat("\n", "-- cluster ", i, " --", "\n", sep="")
# 	cat(as.character(DATA[x[[i]],1]), sep="\n")
# 	cat("\n\n")
# }

# ----------------------------------------
# // Differential Expression //
# see 2018-02-17/compare-across-funcs.R


# ----------------------------------------
