# This script read in interesting gene list file and convert to integer ID
#perl $src_dir/parse_int_multi.pl $interest_file $geneset_use_file\_gene_id
parse_int_multi <- function(interest, geneset_gene_id){
	##geneset_use_file_gene_id input

	##geneset input
	#interest <- read.table(file=interest_file, header = F, sep = "\t", fill = TRUE)
	#interest <- as.matrix(interest)
	
	names(interest) <- c("gene_Name", "sample")
	sample_anno <- unique(interest$sample)
	sample_id <- data.frame(sample=sample_anno, id=1:length(sample_anno))
	tmp1 <- left_join(interest, geneset_gene_id, by = "gene_Name")
	tmp2 <- left_join(tmp1, sample_id, by = "sample")
	
	#colnames
	interest_id <- data.frame(tmp2[,3],tmp2[,4])
	names(interest_id) <- c("gene_id", "sample_id")
	
	diff_gene <- setdiff(interest$gene_Name, geneset_gene_id$gene_Name)
	k = length(diff_gene)

	if(k>0){
		cat("There are", k, "genes from interesting gene list not found in the", "geneset_use_file_gene_id","!\n")
	}

	##output
	out <- list(interest_gene_id=interest_id, interest_gene_sample=sample_id)
	return(out)
}

