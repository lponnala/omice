
# ----------------------------------------
# // Heirarchical Clustering //
# see script.py for preparing the data
# [the following analysis is based on 2014-11-24/analyz.R]

DATA = readr::read_csv("cluster_data.csv")
dim(DATA)
colnames(DATA)

D = DATA[,-1]
sapply(D,class)
unique(sapply(D,class))

# ~~ Dendrogram ~~
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

# ~~ Heatmap ~~
library(gplots)

Dm = as.matrix(D)
# colnames(Dm) = c("WT-FL", "Met1-FL", "WT-NL", "Met1-NL")

# helpful links for heatmap
# https://www.biostars.org/p/14156/
# http://mannheimiagoesprogramming.blogspot.com/2012/06/drawing-heatmaps-in-r-with-heatmap2.html

png(filename = "heatmap.png", width=960, height=480, units="px")
heatmap.2(Dm, Rowv = as.dendrogram(D_hc), dendrogram = "row", col=redgreen(75), scale="row", key=FALSE, density.info="none", trace="none", cexCol=0.8, labRow=NA)
dev.off()


# ----------------------------------------
# // Differential Expression: ANOVA //
# see script.py for preparing the data
# [the following analysis is based on 2018-02-17/compare-across-funcs.R]

y = readr::read_csv("diffexp_anova_data.csv")

y %>% dplyr::count(gen,rep)
y %>% dplyr::count(gen)
y %>% dplyr::count(rep)

# // 2-way anova for testing differences in genotype+function //
y %>% dplyr::group_by(gen,protein_id) %>% dplyr::summarise(m = mean(val)) %>% dplyr::ungroup() %>% dplyr::summarise(cv = sd(m)/mean(m))
# .. do anova ..
ya = y %>% dplyr::mutate(protein_id = as.factor(protein_id), gen = as.factor(gen))
with_interaction = FALSE
if (with_interaction) {
  aov2 = stats::aov(val ~ protein_id*gen, data = ya) # in genotype, function, and their interaction
  summary(aov2) %>% print() # protein_id "effect" is significant
  aov2 %>% broom::tidy() %>% print()
  out_file = "two-way-anova-with-interaction.xlsx"
} else {
  aov2 = stats::aov(val ~ protein_id + gen, data = ya) # in genotype, function
  summary(aov2) %>% print() # protein_id "effect" is very significant
  aov2 %>% broom::tidy() %>% print()
  out_file = "two-way-anova.xlsx"
}
# .. which pairs of proteins are different? ..
diffexp_proteins = stats::TukeyHSD(aov2, which = "protein_id") %>% broom::tidy() %>% tibble::as_tibble() %>% dplyr::arrange(`adj.p.value`)
diffexp_proteins %>% dplyr::filter(`adj.p.value` < 0.01) %>% dplyr::count()

# function to write to excel
toxl = function(data_frames, sheet_names = NA, round = FALSE, column_width = 15, out_file = tempfile(tmpdir = getwd(), fileext = ".xlsx"), open_file = FALSE) {
  # if data_frames is not a list, try to make it one
  if (dplyr::is.tbl(data_frames)) {
    data_frames = list(data_frames)
  } else if (is.data.frame(data_frames)) {
    data_frames = list(as_data_frame(data_frames))
  }
  
  if ((length(sheet_names) == 1) && is.na(sheet_names)) {
    sheet_names = as.character(seq_along(data_frames))
  } else {
    if (length(sheet_names) != length(data_frames)) {
      stop("sheet_names and data_frames are not of same length!")
    }
  }
  
  wb = openxlsx::createWorkbook()
  for (i in seq_along(sheet_names)) {
    openxlsx::addWorksheet(wb, sheetName = sheet_names[[i]])
    openxlsx::writeDataTable(wb, sheet = sheet_names[[i]], x = data_frames[[i]], startCol = 1, startRow = 1)
    openxlsx::setColWidths(wb, sheet = sheet_names[[i]], cols = 1:ncol(data_frames[[i]]), widths = column_width)
    if (round) {
      style = openxlsx::createStyle(numFmt = "#,#", halign = "center", border = "LeftRight")
    } else {
      style = openxlsx::createStyle(halign = "center", border = "LeftRight")
    }
    openxlsx::addStyle(wb, sheet = sheet_names[[i]], style = style, rows = seq(1, nrow(data_frames[[i]]) + 1), cols = 1:ncol(data_frames[[i]]), gridExpand = TRUE)
  }
  openxlsx::saveWorkbook(wb, file = out_file, overwrite = TRUE)
  
  if (open_file) file.show(out_file)
}

# .. write output ..
list(
  aov2 %>% broom::tidy(), 
  stats::TukeyHSD(aov2, which = "protein_id") %>% broom::tidy() %>% dplyr::arrange(`adj.p.value`)
) %>% toxl(sheet_names = c("2-way ANOVA","Compare protein"), open_file = FALSE, out_file = out_file)


# // 1-way anova //
# ~~ testing differences in genotype ~~
y %>% dplyr::group_by(gen) %>% dplyr::summarise(m = mean(val)) %>% dplyr::summarise(cv = sd(m)/mean(m))
aov1 = stats::aov(val ~ gen, data = y)
summary(aov1) # no significant "genotype" effect
aov1 %>% broom::tidy()
# ~~ testing differences in protein_id ~~
y %>% dplyr::group_by(protein_id) %>% dplyr::summarise(m = mean(val)) %>% dplyr::summarise(cv = sd(m)/mean(m))
aov1 = stats::aov(val ~ protein_id, data = y)
summary(aov1) # there is a significant "protein_id" effect
aov1 %>% broom::tidy()


# ----------------------------------------
# // Differential Expression: GLEE //

rm(list = ls(envir = globalenv()), envir = globalenv())

# setup
source("glee-funcs.R") # copied from https://github.com/lponnala/glee-r-pkg/blob/master/gleeR.r
type = c("wt-k1","k1-k6","wt-k6")[3]

D = readr::read_csv(paste0("data_glee_",type,".csv"))
out_stub = paste0("glee_",type)

Data = D
nA = 2
nB = 2
fit_type = "cubic"
num_iter = 10000
num_digits = 4

# run the procedure
Prot = unlist(Data[,1])
A = as.matrix(Data[,1+(1:nA)])
B = as.matrix(Data[,1+nA+(1:nB)])
if (!data_ok(Data,nA,nB)) {
  msg = paste("check input data! spreadsheet must have",
  "(1) the right number of columns",
  "(2) positive finite values",
  "\n",sep="\n")
  stop(msg)
}
m = fit_model(A, B, fit_type)
model_fit_plots(m, outfile=paste0(out_stub,"-fitplots.png"))
stn_pval = calc_stn_pval(A, B, m, num_iter)
stn_pval_plots(stn_pval, outfile=paste0(out_stub,"-stn-pval.png"))
tab = diff_exp_table(stn_pval, Prot, num_digits)

# write out in the same order in which proteins were listed in the input
Data %>% dplyr::select(protein_id) %>%
  dplyr::left_join(tab %>% tibble::as_tibble() %>% dplyr::rename(protein_id = Name), by = "protein_id") %>%
  readr::write_csv(path = paste0(out_stub,"-results.csv"))


