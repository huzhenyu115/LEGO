# This script read in original gene set file and convert to integer ID
# usage: Rscript R/parse_gs.R <original gene set file> <gene Id file>
parse_gs <- function(geneset_file, network_file_id){
	##network_id input
	#network_id
	#network$network_file_id
	
	##geneset input
	#geneset <- read.table(file=geneset_file, header = F, sep = "\t", fill = TRUE)
	#geneset <- as.matrix(geneset)

	geneset <- geneset_file
	network_id <- network_file_id
	max_id <- nrow(network_id)
	
	#geneset_use_file_id output
	path <- unique(geneset[,2])
	gene_id <- data.frame(pathway_Name=path, pathway_Id=1:length(path))
	
	#geneset_use_file_gene_id output
	gene <- unique(geneset[,1])
	diff_gene <- setdiff(gene, network_id$gene_Name)
	geneset_gene_id <- rbind(network_id, data.frame(gene_Name=diff_gene, id=(max_id+1):(max_id+length(diff_gene))))
	
	#geneset_use_file_id_gs output
	annotation1 <- geneset_gene_id
	annotation2 <- gene_id
	names(geneset) <- c("gene_Name", "pathway_Name")
	tmp1 <- left_join(geneset, annotation1, by = "gene_Name")
	tmp2 <- left_join(tmp1, annotation2, by = "pathway_Name")
	#colnames
	gene_pathway_id <- data.frame(tmp2[,3],tmp2[,4])
	names(gene_pathway_id) <- c("gene_id","pathway_id")

	##Number of genes not found in the networ
	#geneset_use_file_gene_id <- paste(geneset_use_file, "gene_id", sep="_")
	k <- nrow(gene_id)-max_id
	if(k>0){
		cat("There are", k, "genes not found in the network;Add them in", "geneset_use_file_gene_id", ";\n")
	}

	##output

	out <- list(geneset_use_file_id=gene_id, geneset_use_file_gene_id=geneset_gene_id, geneset_use_file_id_gs=gene_pathway_id)
	return(out)
}



