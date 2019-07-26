library(xlsx)

R = read.xlsx("mrnaData.xls",sheetIndex=1)
names(R) = c('matched','geneID','protID','cov1','cov4','cov9','cov14','rpkm1','rpkm4','rpkm9','rpkm14')

P = read.xlsx("protData.xls",sheetIndex=1)
names(P) = c('protID','groupID','nsaf1','nsaf4','nsaf9','nsaf14','nadjspc1','nadjspc4','nadjspc9','nadjspc14')

D = merge(R,P,by="protID")
ind = grep(pattern="^GRMZM",D$protID)
D = D[ind,]

save(D, file="matchedGRMZMdataNoGrouping.RData")
write.xlsx(D, file="matchedGRMZMdataNoGrouping.xlsx",row.names=FALSE,showNA=FALSE)


# R = read.xlsx("allMatchedRNA.xlsx",sheetIndex=1)
# names(R) = c('matched','geneID','protID','cov1','cov4','cov9','cov14','rpkm1','rpkm4','rpkm9','rpkm14')
# save(R, file="allMatchedRNA.RData")

# P = read.xlsx("groupIDprot.xlsx",sheetIndex=1)
# names(P) = c('protID','groupID','groupContent','nsaf1','nsaf4','nsaf9','nsaf14','nadjspc1','nadjspc4','nadjspc9','nadjspc14')
# save(P, file="groupIDprot.RData")

# load(file="allMatchedRNA.RData");
# load(file="groupIDprot.RData");

# D = merge(R,P,by="protID")
# save(D, file="matchedDataNoGrouping.RData")
# write.xlsx(D,file="matchedDataNoGrouping.xlsx",row.names=FALSE,showNA=FALSE)

# load(file="matchedDataNoGrouping.RData")
