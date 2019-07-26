
# ==== Compare GLEE and QSPEC results ====

comps = c("prep-wt","oop-wt","triple-wt","prep-triple","oop-triple")

GLEE_pVal_cutoff = 0.01
QSPEC_FDR_cutoff = 0.1

for (comp in comps) {
	cat("\n", "// ", comp, " //", "\n", sep="")
	
	GLEE_results_file = paste0("glee/",comp,"/",comp,".glee-diffexp.csv")
	QSPEC_results_file = paste0("qspec/",comp,"/",comp,".qspec_results.csv")
	
	G = read.csv(GLEE_results_file, header=TRUE)
	Q = read.csv(QSPEC_results_file, header=TRUE)

	G_signif_ind = which(as.numeric(as.character(G$pVal)) < GLEE_pVal_cutoff)
	G_signif_prot = as.character(G[G_signif_ind,1])
	cat("# signif proteins as per GLEE = ", length(G_signif_prot), "\n", sep="")

	Q_signif_ind = which(Q$fdr < QSPEC_FDR_cutoff)
	Q_signif_prot = as.character(Q[Q_signif_ind,1])
	cat("# signif proteins as per QSPEC = ", length(Q_signif_prot), "\n", sep="")

	cat("number of signif proteins identified by either = ", length(union(G_signif_prot, Q_signif_prot)), "\n", sep="")
	cat("number of signif proteins identified by both = ", length(intersect(G_signif_prot, Q_signif_prot)), "\n", sep="")
}

