
infile = "raw-data.txt"
outfile = "datamatrix.txt"
outfile_nolen = "datamatrix-noLength.txt"

y = readLines(infile)
D = vector(length=0)
for (k in 1:8) { D = cbind(D, y[seq((k-1)*11 + 2, k*11)]) }
colnames(D) = y[ (1:8 - 1)*11 + 1 ]

cat(colnames(D), file=outfile, sep="\t", append=F)
cat("\n", file=outfile, append=T)
for (i in 1:dim(D)[1]) {
	cat(D[i,], file=outfile, sep="\t", append=T)
	cat("\n", file=outfile, append=T)
}

D = D[,-2]
cat(colnames(D), file=outfile_nolen, sep="\t", append=F)
cat("\n", file=outfile_nolen, append=T)
for (i in 1:dim(D)[1]) {
	cat(D[i,], file=outfile_nolen, sep="\t", append=T)
	cat("\n", file=outfile_nolen, append=T)
}

# ----------

cat(y[ (1:8 - 1)*11 + 1 ], file=outfile, sep="\t", append=F)
cat("\n", file=outfile, append=T)
for (k in 1:8) {
	cat(y[seq((k-1)*11 + 2, k*11)], file=outfile, sep="\t", append=T)
	cat("\n", file=outfile, append=T)
}


