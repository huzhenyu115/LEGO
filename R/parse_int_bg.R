# This script read in interesting gene list file and convert to integer ID
parse_int_bg <- function(interest, bg_file, geneset_use_file_gene_id, 1){
	##geneset_use_file_gene_id input
	#geneset_use_file_gene_id <- read.table(file=geneset_use_file_gene_id, header = F, sep = "\t", fill = TRUE)
	#geneset_use_file_gene_id <- as.matrix(geneset_use_file_gene_id)
	##bg_file input
	bg <- read.table(bg_file, header = F, sep = "\t", fill = TRUE)
	#bg <- as.matrix(bg)
	
	names(bg) <- "gene_Name"
	tmp1 <- left_join(bg, geneset_gene_id, by = "gene_Name")
	bg_id <- tmp1[,2]
	names(bg_id) <- "gene_id"
	
	#interest_gene_id
	names(interest) <- "gene_Name"
	tmp <- left_join(interest, geneset_gene_id, by = "gene_Name")
	#colnames
	interest_id <- tmp[,2]
	names(interest_id) <- "gene_id"
	
	diff_gene <- setdiff(interest, geneset_gene_id$gene_Name)
	k = length(diff_gene)

	if(k>0){
		cat("There are", k, "genes from interesting gene list not found in the", "geneset_use_file_gene_id","!\n")
	}

	##output	
	out <- list(interest_gene_id=interest_id, bg_gene_id=bg_id)
	return(out)

}
