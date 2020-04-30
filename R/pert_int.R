# label shuffling
# Rscript R/pert_int.r gene_id each perm_times tmp_input
pert_int <- function(gene_id, i, perm_times){
	##input
	nw_file <- gene_id
	size <- i  #151
	per_num <- perm_times  #1000
	# 7SK CLCF1   1
	# A26C1B  1 
	nw_file_input <- read.table(nw_file, header = F, sep = "\t", fill = TRUE)
	#nw_file_input <- as.matrix(nw_file_input)
	all_w <- data.frame(sample(nw_file_input[,1]))
	##
	all_result <- data.frame()
	for(j in 1:per_num){
		all_w <- data.frame(sample(all_w))
		tmp <- all_w[1:size,]
		all_tmp <- data.frame(gene_Name=tmp, number=rep(j,length(tmp)))
		all_result <- rbind(all_result, all_tmp)
	}
	names(all_result) <- c("gene_Name", "number")
	#write.table(all_result, file=tmp_input, sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
	
	out <- all_result
	return(out)
}

