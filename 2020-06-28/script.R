
# # --------------------------------------------------------------------------------
# # do PCA
# # --------------------------------------------------------------------------------
# rm(list = ls(envir = globalenv()), envir = globalenv())
# load("N.RData")
# method = c("prcomp","princomp")[[1]]
# 
# if (method == "prcomp") {
#   y = prcomp(N[,-1])
#   pc = y$rotation
# } else if (method == "princomp") {
#   y = princomp(N[,-1])
#   pc = y$loadings
# }
# 
# print(summary(y))
# print(pc)
# # print(pc, cutoff=0.0001)
# 
# jpeg('NadjSPC-pcaPlot.jpg')
# plot(pc[1:3,1:2],type='p',pch=15,main="PCA using NadjSPC data",xlab='PC1',ylab='PC2',xlim=c(min(pc[,1])-0.1,max(pc[,1])+0.1),ylim=c(min(pc[,2])-0.1,max(pc[,2])+0.1))
# points(pc[4:6,1:2],type='p',pch=17)
# legend('topleft',inset=0.05,legend=c('wt','cgep'),pch=c(15,17),col='black')
# dev.off()
# 
# write.csv(pc[,1:2], file = "NadjSPC-pca.csv")
# # --------------------------------------------------------------------------------


# # --------------------------------------------------------------------------------
# # do t-test
# # --------------------------------------------------------------------------------
# rm(list = ls(envir = globalenv()), envir = globalenv())
# dato = readxl::read_excel("t-test-forLalit.xlsx")
# 
# dat = dato %>% dplyr::select(contains("func"),contains("nadj"))
# colnames(dat) %>% print()
# colnames(dat) = c("simple functional name","wt rep1","wt rep2","wt rep3","cgep rep1","cgep rep2","cgep rep3")
# 
# dat %<>% 
#   tidyr::gather(key = "key", value = "val", -`simple functional name`) %>% 
#   tidyr::separate(key, into = c("gen","rep"), sep = " ") %>% 
#   dplyr::group_by(`simple functional name`) %>% tidyr::nest() %>% 
#   dplyr::mutate(ttN = purrr::map(data, function(x) t.test(val ~ gen, data = x, var.equal = FALSE, paired = FALSE, alternative = "two.sided"))) %>% 
#   dplyr::mutate(ttNP = purrr::map(ttN, function(x) x %>% broom::tidy() %>% dplyr::transmute(statistic,`p.value`,inference = ifelse(`p.value` < 0.05, "NOTEQ", "EQUAL"),method,alternative))) %>% 
#   tidyr::unnest(ttNP) %>% 
#   dplyr::select_if(~ !is.list(.x))
# 
# dato %>% dplyr::left_join(dat, by = "simple functional name") %>% 
#   readr::write_csv("ttest-results.csv")
# 
# # inspect the results (copy-paste the t-test columns into Klaas's original spreadsheet)
# file.show("ttest-results.csv")
# # --------------------------------------------------------------------------------


# # --------------------------------------------------------------------------------
# # compare QSPEC vs GLEE output
# # --------------------------------------------------------------------------------
# rm(list = ls(envir = globalenv()), envir = globalenv())
# glee = readr::read_csv("glee-results.csv")
# qspec = readr::read_csv("qspec-results.csv")
# 
# symdiff(qspec$Protein,qspec$Protein)
# 
# glee_cutoff = 0.01
# qspec_cutoff = 0.05
# 
# glee %>% dplyr::filter(pVal < glee_cutoff) %>% nrow()
# qspec %>% dplyr::filter(fdr < qspec_cutoff) %>% nrow()
# 
# intersect(
#   glee %>% dplyr::filter(pVal < glee_cutoff) %>% .$Accession,
#   qspec %>% dplyr::filter(fdr < qspec_cutoff) %>% .$Protein
# ) %>% length()
# # --------------------------------------------------------------------------------


# # --------------------------------------------------------------------------------
# # run QSPEC
# # --------------------------------------------------------------------------------
# rm(list = ls(envir = globalenv()), envir = globalenv())
# load("A.RData")
# datfile = "datamatrix.txt"
# outfile = "qspec-results.csv"
# 
# # write the data in proper format
# cat("Protein\tLength\t0\t0\t0\t1\t1\t1\n", file=datfile, sep="", append=FALSE)
# A %>% dplyr::mutate(Length = 1) %>%
#   dplyr::select(Accession, Length, everything()) %>%
#   readr::write_tsv(path = datfile, append = TRUE)
# # inspect the file
# file.show(datfile)
# 
# # run QSPEC on Linode
# system("rm -f *qspec*")
# system(paste0("../qprot-linux/qprot_1.3.0/qspec-param ",datfile," 5000 20000 0"))
# system(paste0("../qprot-linux/qprot_1.3.0/getfdr ",datfile,"_qspec"))
# 
# # clean up the output (remove the "dummy" length column, and assign proper column names)
# R = readr::read_tsv(file = paste0(datfile,"_qspec_fdr"), col_names = TRUE)
# R %<>% dplyr::select(-Len)
# colnames(R)[2:7] = colnames(A)[2:7]
# R %>% readr::write_csv(path = outfile)
# 
# # inspect the output (copy-paste it into Klaas's original spreadsheet, ensure proteins are in same order via spot check)
# file.show(outfile)
# # --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# run GLEE
# --------------------------------------------------------------------------------
rm(list = ls(envir = globalenv()), envir = globalenv())

# setup
source("glee-funcs.R") # copied from https://github.com/lponnala/glee-r-pkg/blob/master/gleeR.r
load("N.RData")
Data = N
nA = 3
nB = 3
fit_type = "cubic"
num_iter = 10000
num_digits = 4
out_stub = "glee-NadjSPC" # file = "glee-results-NadjSPC.csv"

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
Data %>% dplyr::select(Accession) %>%
  dplyr::left_join(tab %>% tibble::as_tibble() %>% dplyr::rename(Accession = Name), by = "Accession") %>%
  readr::write_csv(path = paste0(out_stub,"-results.csv"))

# inspect the output (copy-paste it into Klaas's original spreadsheet, ensure proteins are in same order via spot check)
file.show(paste0(out_stub,"-results.csv"))
# --------------------------------------------------------------------------------


# # --------------------------------------
# # run correlations
# # --------------------------------------
# rm(list = ls(envir = globalenv()), envir = globalenv())
# load("A.RData")
# cor(A[,-1]) %>% print()
# # write.csv(cor(A[,-1]), file="corr-adjSPC.csv")
# load("N.RData")
# cor(N[,-1]) %>% print()
# # write.csv(cor(N[,-1]), file="corr-NadjSPC.csv")
# # --------------------------------------


# # --------------------------------------------------------------------------------
# # prepare data
# # --------------------------------------------------------------------------------
# rm(list = ls(envir = globalenv()), envir = globalenv())
# 
# dat = readxl::read_excel("forLalit-Reiset-June2020v2.xlsx")
# colnames(dat) %>% cat(sep="\n")
# 
# A = dat %>% dplyr::select(c(1,3:8))
# colnames(A) %>% cat(sep="\n")
# colnames(A) = c("Accession","wt_rep1","wt_rep2","wt_rep3","trap_rep1","trap_rep2","trap_rep3")
# dim(A)
# A %<>% dplyr::filter(!is.na(Accession)) # nested "missing" rows in spreadsheet
# dim(A)
# A %>% dplyr::filter(is.finite(rowSums(.[-1]))) %>% dim()
# A %>% dplyr::mutate_if(is.numeric, ~ ifelse(is.finite(.x),.x,0)) %>% dplyr::filter(is.finite(rowSums(.[-1]))) %>% dim()
# A %<>% dplyr::mutate_if(is.numeric, ~ ifelse(is.finite(.x),.x,0))
# A %>% dplyr::filter(is.finite(rowSums(.[-1]))) %>% dim()
# A %<>% dplyr::mutate(Accession = stringr::str_replace_all(Accession," ","_")) # QSPEC does not like spaced in protein-name
# save(A, file = "A.RData")
# 
# N = dat %>% dplyr::select(c(1,10:15))
# colnames(N) %>% cat(sep="\n")
# colnames(N) = c("Accession","wt_rep1","wt_rep2","wt_rep3","trap_rep1","trap_rep2","trap_rep3")
# N %<>% dplyr::filter(!is.na(Accession))
# dim(N)
# N %<>% dplyr::mutate(Accession = stringr::str_replace_all(Accession," ","_")) # to be consistent with what we did for adjSPC above
# save(N, file = "N.RData")
# 
# # // normalized to ClpC1 //
# 
# A = dat %>% dplyr::select(c(1,17:22))
# colnames(A) %>% cat(sep="\n")
# colnames(A) = c("Accession","wt_clpc1_rep1","wt_clpc1_rep2","wt_clpc1_rep3","trap_clpc1_rep1","trap_clpc1_rep2","trap_clpc1_rep3")
# dim(A)
# A %<>% dplyr::filter(!is.na(Accession)) # nested "missing" rows in spreadsheet
# dim(A)
# A %>% dplyr::filter(is.finite(rowSums(.[-1]))) %>% dim()
# A %>% dplyr::mutate_if(is.numeric, ~ ifelse(is.finite(.x),.x,0)) %>% dplyr::filter(is.finite(rowSums(.[-1]))) %>% dim()
# A %<>% dplyr::mutate_if(is.numeric, ~ ifelse(is.finite(.x),.x,0))
# A %>% dplyr::filter(is.finite(rowSums(.[-1]))) %>% dim()
# A %<>% dplyr::mutate(Accession = stringr::str_replace_all(Accession," ","_")) # QSPEC does not like spaced in protein-name
# save(A, file = "A_clpc1.RData")
# 
# N = dat %>% dplyr::select(c(1,24:29))
# colnames(N) %>% cat(sep="\n")
# colnames(N) = c("Accession","wt_clpc1_rep1","wt_clpc1_rep2","wt_clpc1_rep3","trap_clpc1_rep1","trap_clpc1_rep2","trap_clpc1_rep3")
# N %<>% dplyr::filter(!is.na(Accession))
# dim(N)
# N %<>% dplyr::mutate(Accession = stringr::str_replace_all(Accession," ","_")) # to be consistent with what we did for adjSPC above
# save(N, file = "N_clpc1.RData")
# # --------------------------------------------------------------------------------
