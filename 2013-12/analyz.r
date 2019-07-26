
# R = read.xlsx("mrna.xlsx",sheetIndex=1)

# P = read.xlsx("prot.xlsx",sheetIndex=1)
# names(P)[1]="protID"

# D = merge(R,P,by="protID")
# gid = which(is.finite(D$groupID))
# nogid = which(!is.finite(D$groupID))

# save(D,file="allData.RData")
# write.xlsx(D,file="allData.xlsx",row.names=FALSE,showNA=FALSE)
# write.xlsx(D[gid,],file="gidData.xlsx",row.names=FALSE,showNA=FALSE)
# write.xlsx(D[nogid,],file="nogidData.xlsx",row.names=FALSE,showNA=FALSE)

load(file="allData.RData")

covInd = 4:7
rpkmInd = 8:11
nsafInd = 13:16
nadjspcInd = 17:20

# using average across sections
cat('RPKM & NSAF\tCorrelation\tPvalue\n')
y=cor.test(log10(rowMeans(D[,rpkmInd])),log10(rowMeans(D[,nsafInd])),method="pearson")
cat('Pearson (Rp)',round(y$estimate,6),round(y$p.value,6),'\n',sep='\t')
y=cor.test(log10(rowMeans(D[,rpkmInd])),log10(rowMeans(D[,nsafInd])),method="spearman")
cat('Spearman (Rs)',round(y$estimate,6),round(y$p.value,6),'\n',sep='\t')

cat('COV & NadjSPC\tCorrelation\tPvalue\n')
y=cor.test(log10(rowMeans(D[,covInd])),log10(rowMeans(D[,nadjspcInd])),method="pearson")
cat('Pearson (Rp)',round(y$estimate,6),round(y$p.value,6),'\n',sep='\t')
y=cor.test(log10(rowMeans(D[,covInd])),log10(rowMeans(D[,nadjspcInd])),method="spearman")
cat('Spearman (Rs)',round(y$estimate,6),round(y$p.value,6),'\n',sep='\t')



