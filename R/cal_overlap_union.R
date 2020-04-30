cal_overlap_union <- function(geneset_use_file){
	# BIOCARTA_RELA_PATHWAY   BIOCARTA_AKT_PATHWAY    Gene_set_overlap    0.37
	##geneset_use_file input
	#geneset_use_file

	all_gs <- unique(geneset_use_file$pathway_Name)
	gs_num <- length(all_gs)
	total_adj <- gs_num*(gs_num-1)/2
	cat(total_adj, "\n")

	result <- data.frame()
	for(i in 1:(gs_num-1)){
		gs1 <- all_gs[i]
		all_gene1 <- filter(geneset_use_file, pathway_Name==gs1)[1]
		B <- nrow(all_gene1)
		for(j in (i+1):gs_num){
			gs2 <- all_gs[j]
			all_gene2 <- filter(geneset_use_file, pathway_Name==gs2)[1]
			C <- nrow(all_gene2)
			D <- nrow(intersect(all_gene1,all_gene2))
			if(D == 0){next}
			A <- B+C-D
			jac <- D/A
			result <- rbind(result, data.frame(gs1, gs2, jac))
		}
	}
	##output
	write.table(result, file=paste("geneset_file", "overlap_union.txt", sep="_"), sep = "\t", quote = FALSE, col.names = FALSE)
	out <- result
	return(out)

}
