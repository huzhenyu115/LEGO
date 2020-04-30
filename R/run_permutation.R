# run_permutation(interest_file_id.out, geneset_use_file, network_file, geneset_file, main_dir, perm_times)
# Rscript R/run_permutation.r interest_file_id.out /mid_data_dir multi
run_permutation <- function(interest_file_id.out
							,geneset_use_file_name
							,network_file
							,geneset_file
							,main_dir
							,perm_times)
{
	input <- interest_file_id.out
	geneset <- geneset_use_file_name
	net_file <- network_file
	geneset_file <- geneset_file
	main_dir <- main_dir
	perm_times <- perm_times  #1000
	gene_id <- paste(net_file, "id", sep="_")
	mid_data_dir <- paste(geneset, "_mid_data")
	## check dir
	if(file.exists(mid_data_dir)){
		next
	}else{
		dir.create(mid_data_dir)
	}
	##input
	interest_file_id_out <- read.table(file=input, header = T, sep = "\t", fill = TRUE)
	#interest_file_id_out <- as.matrix(interest_file_id_out)
	each <- sort(unique(interest_file_id_out$intSize_net))
	
	for(i in each){	
		final_file <- paste(mid_data_dir, "/", i, ".RData", sep="")
		if(file.exists(final_file)){
			next
		}
		
		cat("Rscript pert_int.R", "gene_id", i, perm_times, "tmp_input", "\n")
		tmp_input <- pert_int(gene_id, i, perm_times)
		cat("Rscript LEGO_noperm.R", net_file, geneset_file, tmp_input, "-multi 1 -noR 0 -min 0 -max 100000 \n")
		LEGO_noperm(net_file, geneset_file, tmp_input, multi=1, noR=0, min=0, max=100000)
		tmp_input_id.out <- paste(tmp_input, "id.out", sep="_")
		cat("Rscript summary_pert.R", tmp_input_id.out, final_file, main_dir, "\n")
		summary_pert(tmp_input_id.out, final_file, main_dir)
		lf<-list.files(main_dir,pattern = tmp_input)
		file.remove(lf) 
	}
	cat("Input size:", i, "\n")
	
	out <- 
	return()
}

