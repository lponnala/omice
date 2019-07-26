library(xlsx)

XLSXFILE = "clusterset-forLalit-sept162014.xlsx"
DATAFILE = "clusterset-forLalit-sept162014.RData"

# colNames = unlist(read.xlsx(XLSXFILE, sheetIndex=1, rowIndex=1, colIndex=1:26, as.data.frame=FALSE, header=FALSE))
# DATA = read.xlsx(XLSXFILE, sheetIndex=1, rowIndex=2:3520, colIndex=1:26, as.data.frame=TRUE, header=FALSE)
# names(DATA) = colNames
# save(DATA, file=DATAFILE)

load(DATAFILE)
thresh = 1e-4
save_plots = TRUE

W = DATA[DATA[,15] >= thresh, 17:21]
W_dd = as.dist((1-cor(t(W)))/2)
W_hc = hclust(W_dd, method="average")
if (save_plots) { png(filename = "Wt.png", width=960, height=480, units="px") }
plot(W_hc, labels=FALSE, hang=-1, frame.plot=FALSE, main="Wild-type", sub="", xlab="", ylab="correlation-based distance")
num_clust = 6
rect.hclust(W_hc, k=num_clust, border = 1 + 1:num_clust)
if (save_plots) {	dev.off() }

M = DATA[DATA[,16] >= thresh, 22:26]
M_dd = as.dist((1-cor(t(M)))/2)
M_hc = hclust(M_dd, method="average")
if (save_plots) { png(filename = "M.png", width=960, height=480, units="px") }
plot(M_hc, labels=FALSE, hang=-1, frame.plot=FALSE, main="Mutant", sub="", xlab="", ylab="correlation-based distance")
num_clust = 7
rect.hclust(M_hc, k=num_clust, border = 1 + 1:num_clust)
if (save_plots) {	dev.off() }

WM = DATA[DATA[,15] >= thresh & DATA[,16] >= thresh, 17:26]
WM_dd = as.dist((1-cor(t(WM)))/2)
WM_hc = hclust(WM_dd, method="average")
if (save_plots) { png(filename = "Wt_M_both.png", width=960, height=480, units="px") }
plot(WM_hc, labels=FALSE, hang=-1, frame.plot=FALSE, main="Wild-type and Mutant (threshold on both)", sub="", xlab="", ylab="correlation-based distance")
num_clust = 9
rect.hclust(WM_hc, k=num_clust, border = 1 + 1:num_clust)
if (save_plots) {	dev.off() }

WMa = DATA[DATA[,14] >= thresh, 17:26]
WMa_dd = as.dist((1-cor(t(WMa)))/2)
WMa_hc = hclust(WMa_dd, method="average")
if (save_plots) { png(filename = "Wt_M_avg.png", width=960, height=480, units="px") }
plot(WMa_hc, labels=FALSE, hang=-1, frame.plot=FALSE, main="Wild-type and Mutant (threshold on average)", sub="", xlab="", ylab="correlation-based distance")
num_clust = 9
rect.hclust(WMa_hc, k=num_clust, border = 1 + 1:num_clust)
if (save_plots) {	dev.off() }


# # -- to identify clusters manually --
# x = identify(hc)
# DATA[x[[1]],1]

