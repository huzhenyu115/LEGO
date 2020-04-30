############################################## LEGO
# Rscript LEGO.R <network file> <geneset file> <interest file> [options]
# -pre: <whether or not to pre-run to generate mid files,1-yes,0-no(e.g:0)>
# -multi: <multi or not:1-yes,0-no(e.g:0)>
# -bgNE: <bgNE for the network,e.g:0.25>
# -min: <minimun gene set size,e.g:5>
# -max: <maximum gene set size,e.g:10000>
# -adj: <adjusted methods,e.g:fdr>
# -p_thre: <p value cutoff,e.g:0.05>
# -fisher: <whether or not run fisher exact test, 1-yes,0-no,e.g: 0>
# -filter: <whether or not run result cluster and filter step, 1-yes, 0-no, e.g: 0>
# -bg_file: <background file, e.g: bg.txt, if do not have background file, leave it blank> \n";
# -perm_times: <permutation time,default: 1000>

#network_file
#geneset_file
#interest_file
#network_file_id
#network_file_id_net
#

LEGO <- function(network_file  ## network file
				         ,geneset_file ## gene set file
                 ,interest_file ## interesting file
                 ,pre_run = 0
                 ,min=5
                 ,max=10000
                 ,multi = 0
                 ,bgNE = 0.25
                 ,adj_meth = "fdr"
                 ,p_thre = 0.05
                 ,fisher = 0
                 ,filter = 0
                 ,bg_file = ""
                 ,perm_times = 1000)
{
	if (missing(network_file)) {
		stop('Please provide "network_file"!')
	}
	if (missing(geneset_file)) {
		stop('Please provide "geneset_file"!')
	}
  if (missing(interest_file)) {
		stop('Please provide "interest_file"!')
  }
	# run program
	# install.packages("Rcpp")
	library(Rcpp)
	library(dplyr)
	cat("\n###### Run begin:", date(), "######\n")
	geneset_use_file <- geneset_file
	interest_file <- read.table(file=interest_file, header=F, sep="\t")
	#interest_file <- as.matrix(interest_file)
  load("data/GeneSet_human.txt_FC2_human_id_gs.rda")
	cat("Run pre-programs to generate mid-files:\n")
	if(pre_run){
		network <- parse_net(network_file)
		#network$network_file_id
		#network$network_file_id_net

		geneset <- parse_gs(geneset_file, network$network_file_id)
		#geneset$geneset_use_file_id
		#geneset$geneset_use_file_gene_id
		#geneset$geneset_use_file_id_gs

		cat("Run LEGO_pre to generate mid files:\n\t", 'LEGO_pre', "network_file_id_net", "geneset_use_file_id_gs", bgNE, "\n")
		system("LEGO_pre network_file_id_net geneset_use_file_id_gs bgNE")
	}else{
		if(!exists("network")){
			network <- parse_net(network_file)
			geneset <- parse_gs(geneset_file, network$network_file_id)

			cat("Run LEGO_pre to generate mid files:\n\t", 'LEGO_pre', "network_file_id_net", "geneset_use_file_id_gs", bgNE, "\n")
			system("LEGO_pre network_file_id_net geneset_use_file_id_gs bgNE")
		}else{
			if(!exists("geneset_use_file_id_gs")){
				geneset <- parse_gs(geneset_file, network$network_file_id)

				system("LEGO_pre network_file_id_net geneset_use_file_id_gs bgNE")
			}
			if(!exists("geneset_use_file_id_gs_NW")){
				system("LEGO_pre network_file_id_net geneset_use_file_id_gs bgNE")
			}
		}
	}
	## no bg
	if(multi){
		interest_gene <- parse_int_multi(interest_file, geneset$geneset_use_file_gene_id)
		#interest_gene$interest_gene_id
		#interest_gene$interest_gene_sample

		cat("Run main program:")
		system("LEGO_mul geneset_use_file_id_gs geneset_use_file_id_gs_NW geneset_use_file_id_gs_CS geneset_use_file_id_gs_GS interest_gene$interest_gene_id min max")
	}else{
		interest_gene_id <- parse_int(interest_file, geneset$geneset_use_file_gene_id)
		cmd = "LEGO geneset_use_file_id_gs geneset_use_file_id_gs_NW geneset_use_file_id_gs_CS geneset_use_file_id_gs_GS interest_gene_id min max"
		cat("\n*Run main program:\n\t", cmd, "\n")
		system(cmd)
		#interest_file_id.out
	}
	## check permutation
	mid_prefix = ""
	bg = c()
	if(bg_file == ""){
		mid_prefix = "mid_data/"
		cat("\n*Check permutation:\n\trun_permutation.R", "interest_file_id.out", "geneset_use_file", "network_file", "geneset_file", perm_times, "\n")
		interest_file_id.out <- run_permutation(interest_file_id.out, geneset_use_file, network_file, geneset_file, main_dir, perm_times)
	}else{
		bg = strsplit(bg_file, "/")
		mid_prefix = paste("mid_data/", bg[[1]][2], "_", sep= "_")
		cat("\n*Check permutation:\n\trun_permutation_bg.R", "interest_file_id.out", "geneset_use_file", "network_file", "geneset_file", main_dir, bg_file, perm_times, "\n")
		interest_file_id.out <- run_permutation_bg(interest_file_id.out, geneset_use_file, network_file, geneset_file, main_dir, bg_file, perm_times)
	}
	##
	if(fisher){
		cat("\n*Run Fisher's Exact test:\n\tRscript enrich.R", "geneset_use_file", "interest_file", "interest_file_fisher", max, min, 1, p_thre, adj_meth, multi, "\n")
		interest_file_fisher <- enrich(database_file=geneset_use_file, input_file=interest_file)
	}
	cat("Rscript interpret_gpd.R", "interest_file_id", "geneset_use_file", adj_meth, p_thre, multi, mid_prefix, "\n")
	input_LEGO <- interpret_gpd(interest_file_id,  geneset_use_file, adj_meth, p_thre, multi, mid_prefix)
	### filter
	if(filter){
		if(!exists("geneset_use_file_overlap_union")){
			cat("\n*Run cal_overlap to calculate overlap value between gene set pairs:\n\tRscript cal_overlap_union.R\n")
		  geneset_use_file_overlap_union <- cal_overlap_union(geneset_use_file)

		}

		cat("\n*Run ORA_filter to get filtered and clustered enriched gene sets:\n\tRscript ORA_filter.R\n")
		interest_file_id.out_LEGO_filter <- ORA_filter(interest_file_id.out_LEGO.txt, interest_file_id.out_LEGO_filter, geneset_use_file, multi)
		if(fisher){
			interest_file_fisher_filter <- ORA_filter(interest_file_fisher, interest_file_fisher_filter, geneset_use_file, multi, 7)
		}
	}
	### finish
	cat("###### Run end:", date(), "######\n\n")
}
