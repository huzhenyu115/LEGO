# This script read in interesting gene list file and convert to integer ID
# $src_dir/parse_int_multi_bg.pl $interest_file $bg_file $geneset_use_file\_gene_id 1
parse_int_multi_bg <- function(interest_file, bg_file, geneset_use_file_gene_id, 1){
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
	out <- list(interest_gene_id=interest_id, interest_gene_sample=sample_id, bg_gene_id=bg_id)
	return(out)

}

