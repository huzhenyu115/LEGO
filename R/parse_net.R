# usage: Rscript R/parse_net.r <original network file>
parse_net <- function(network_file){

	##input
	##network_file
	##FC2_human.txt_id
	x <- network_file
	col_merge <- union(x[,1], x[,2])
	annotation <- data.frame(gene_Name=col_merge, id=1:length(col_merge))
	annotation1 <- data.frame(gene_Name_1=col_merge, id=1:length(col_merge))
	annotation2 <- data.frame(gene_Name_2=col_merge, id=1:length(col_merge))
	tmp1 <- left_join(x, annotation1, by = "gene_Name_1")
	tmp2 <- left_join(tmp1, annotation2, by = "gene_Name_2")
	#FC2_human.txt_id_net
	id_net <- data.frame(tmp2[,4],tmp2[,5],tmp2[,2])
	names(id_net) <- c("id_1","id_2","weight")
	
	tmp3 <- data.frame(x[,1],x[,2])
	del <- nrow(tmp3)-nrow(distinct(tmp3))
	cat("Delete duplicate edges:", del, "\n")
	
	out <- list(network_file_id=annotation, network_file_id_net=id_net)
	
	return(out)
}
