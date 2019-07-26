library(xlsx)

# ==== SECOND PASS ====

XLSFILE = "forlalit-aug15-groupsum.xlsx"
SAVFILE = "forlalit-aug15-groupsum.RData"
columns = 1:26
rows = c(1,1662)
OUTFILE = paste("output.",XLSFILE,sep='')

# heads = unlist(read.xlsx(XLSFILE, sheetIndex=1, startRow=1, colIndex=columns, endRow=1, as.data.frame=FALSE, header=FALSE))
# names(heads) = NULL
# D = read.xlsx(XLSFILE, sheetIndex=1, startRow=rows[1], colIndex=columns, endRow=rows[2], as.data.frame=TRUE, header=TRUE)
# names(D) = heads
# save(D, file=SAVFILE)

load(SAVFILE)

S = data.frame()
for (k in sort(unique(D[,2]))) {
	rows = which(D[,2]==k)
	G = D[rows,]
	
	I = which(G[,14]==max(G[,14]))
	if (length(I)==1) {
		S = rbind(S, cbind(G[I,1:13],as.data.frame(t(colSums(G[,14:26])))))
	} else {
		# pick the alphabetically sorted first accession
		Gmax = G[I,]
		ind = order(Gmax[,1])
		S = rbind(S, cbind(Gmax[ind[1],1:13],as.data.frame(t(colSums(G[,14:26])))))
	}
}

if (file.exists(OUTFILE)) { unlink(OUTFILE) }
write.xlsx(S, file=OUTFILE, col.names=TRUE, row.names=FALSE, append=FALSE, showNA=FALSE)


# # ==== FIRST PASS ====

# XLSFILE = "forlalit-sum-groups&annotate.xlsx"
# SAVFILE = "forlalit-sum-groups&annotate.RData"
# columns = c(1:23, 25:34, 36:45)
# OUTFILE = "output.forlalit-sum-groups&annotate.xlsx"

# # heads = unlist(read.xlsx(XLSFILE, sheetIndex=1, startRow=1, colIndex=columns, endRow=1, as.data.frame=FALSE, header=FALSE))
# # names(heads) = NULL
# # D = read.xlsx(XLSFILE, sheetIndex=1, startRow=1, colIndex=columns, endRow=2001, as.data.frame=TRUE, header=TRUE)
# # names(D) = heads
# # save(D, file=SAVFILE)

# load(SAVFILE)

# S = data.frame()
# for (k in sort(unique(D[,2]))) {
	# rows = which(D[,2]==k)
	# if (length(rows)==0) { stop('cannot identify rows for group ',k,'\n') }
	# G = D[rows,]
	
	# I = which(G[,11]==max(G[,11]))
	# if (length(I)==1) {
		# S = rbind(S, cbind(G[I,1:10],as.data.frame(t(colSums(G[,11:43])))))
	# } else {
		# # pick the alphabetically sorted first accession
		# best = G[I,]
		# ind = order(best[,1])
		# S = rbind(S, cbind(best[ind[1],1:10],as.data.frame(t(colSums(G[,11:43])))))
	# }
# }

# if (file.exists(OUTFILE)) { unlink(OUTFILE) }
# write.xlsx(S, file=OUTFILE, col.names=TRUE, row.names=FALSE, append=FALSE, showNA=FALSE)


# # ==== OLD CODE ====

# XLSFILE = "leafgradientmutant-wtaug2014.xlsx"
# SAVFILE = "leafgradientmutant-wtaug2014.RData"
# columns = c(1:30,32:46)
# OUTFILE = "output.xlsx"

# # heads = unlist(read.xlsx(XLSFILE, sheetIndex=1, startRow=1, colIndex=columns, endRow=1, as.data.frame=FALSE, header=FALSE))
# # D = read.xlsx(XLSFILE, sheetIndex=1, startRow=1, colIndex=columns, endRow=2001, as.data.frame=TRUE, header=TRUE)
# # names(D) = heads
# # save(D, file=SAVFILE)

# load(SAVFILE)

# S = data.frame()
# for (k in unique(D[,2])) {
	# rows = which(D[,2]==k)
	# if (length(rows)==0) { stop('cannot identify rows for group ',k,'\n') }
	# G = D[rows,]
	# best_ind = which.max(G[,13])
	# S = rbind(S, cbind(G[best_ind,1:12],as.data.frame(t(colSums(G[,13:45])))))
# }

# if (file.exists(OUTFILE)) { unlink(OUTFILE) }
# write.xlsx(S, file=OUTFILE, col.names=TRUE, row.names=FALSE, append=FALSE, showNA=FALSE)
