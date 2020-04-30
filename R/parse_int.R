# This script read in interesting gene list file and convert to integer ID
parse_int <- function(interest, geneset_gene_id){
	##geneset_use_file_gene_id input

	##geneset input
	#interest <- read.table(file=interest_file, header = F, sep = "\t", fill = TRUE)
	#interest <- as.matrix(interest)

	names(interest) <- "gene_Name"
	tmp <- left_join(interest, geneset_gene_id, by = "gene_Name")
	#colnames
	interest_id <- tmp[2]
	names(interest_id) <- "gene_id"

	diff_gene <- setdiff(interest$gene_Name, geneset_gene_id$gene_Name)
	k = length(diff_gene)

	if(k>0){
		cat("There are", k, "genes from interesting gene list not found in the", "geneset_use_file_gene_id","!\n")
	}

	##output
	#interest_file_id <- paste(interest_file, "_id")
	#write.table(all_keys3, file=interest_file_id, sep = "\t", quote = FALSE, row.names= FALSE, col.names = FALSE)

	out <- interest_id
	return(out)
}

