############################################## LEGO_noperm
# perl LEGO.pl <network file> <geneset file> <interest file> [options]
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
LEGO_noperm <- function(network_file ## network file
				 ,geneset_file ## gene set file
                 ,interest_file ## interesting file
                 ,main_dir
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
                 ,perm_times = 1000
				 ,no_R = 1)
{
	# create the root directory
	if( !file.exists(main_dir))
    {
		if( !dir.create(main_dir, recursive = TRUE) )
		{
			stop("Unable to create root directory.\n")
		}
    }
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
	cat("\n###### Run begin:", date(), "######\n")
	#network_mark = strsplit(network_file, ".txt")
	#geneset_use_file = paste(geneset_file, network_mark, sep= "_")

	network_file <- read.table(file=network_file, header=F, sep="\t")
	#network_file <- as.matrix(network_file)
	names(network_file)<-c("gene_Name_1","gene_Name_2","weight")
	geneset_file <- read.table(file=geneset_file, header=F, sep="\t")
	#geneset_file <- as.matrix(geneset_file)
	geneset_use_file <- geneset_file
	interest_file <- read.table(file=interest_file, header=F, sep="\t")
	#interest_file <- as.matrix(interest_file)

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
	if(bg_file){
		if(multi){
			interest_gene <- parse_int_multi_bg(interest_file, bg_file, geneset$geneset_use_file_gene_id, 1)
			#interest_gene$interest_gene_id
			#interest_gene$interest_gene_sample
			#interest_gene$bg_gene_id

			cmd = "LEGO_mul_bg geneset_use_file_id_gs geneset_use_file_id_gs_NW geneset_use_file_id_gs_CS geneset_use_file_id_gs_GS interest_file_id min max bg_file_id")
			cat("\n*Run main program:\n\t", cmd, "\n")
			system(cmd)
		}else{
			interest_gene <- parse_int_bg(interest_file, bg_file, geneset$geneset_use_file_gene_id, 1)
			cmd = "LEGO_bg geneset_use_file_id_gs geneset_use_file_id_gs_NW geneset_use_file_id_gs_CS $geneset_use_file_id_gs_GS interest_file_id min max bg_file_id"
			cat("\n*Run main program:\n\t", cmd, "\n")
			system(cmd)
			#interest_file_id.out
		}
		if(fisher){
			cat("\n*Run Fisher's Exact test:\n\tRscript enrich.R\n")
			enrich(geneset_use_file, interest_file, paste(interest_file, "fisher.txt", sep="_"), max, min, 1, p_thre, adj_meth, multi, bg_file)
		}
	}else{
		if(multi){
			interest_gene <- parse_int_multi(interest_file, geneset$geneset_use_file_gene_id)
			cmd = "LEGO_mul geneset_use_file_id_gs geneset_use_file_id_gs_NW geneset_use_file_id_gs_CS geneset_use_file_id_gs_GS interest_file_id min max"
			cat("\n*Run main program:\n\t", cmd, "\n")
			system(cmd)
		}else{
			interest_gene <- parse_int(interest_file, geneset$geneset_use_file_gene_id)
			cmd = "LEGO geneset_use_file_id_gs geneset_use_file_id_gs_NW geneset_use_file_id_gs_CS geneset_use_file_id_gs_GS interest_file_id min max"
			cat("\n*Run main program:\n\t", cmd, "\n")
			system(cmd)
		}
		if(fisher){
			cat("\n*Run Fisher's Exact test:\n\tRscript enrich.R",geneset_use_file, interest_file, paste(interest_file, "fisher.txt", sep="_"), "-max", max, "-min", min, "-list 1 -thre", p_thre, "-adj", adj_meth, "-multi", multi, "\n")
			interest_file_fisher <- enrich(geneset_use_file, interest_file, paste(interest_file, "fisher.txt", sep="_"), max, min, 1, p_thre, adj_meth, multi)
		}
	}
	if(no_R==0){
		break
	}
	input_LEGO <- interpret(interest_file_id.out, geneset_use_file_id, adj_meth, p_thre, multi)
	### filter
	if(filter){
		if(!exists("geneset_use_file_overlap_union")){
			cat("\n*Run cal_overlap to calculate overlap value between gene set pairs:\n\tRscript cal_overlap_noperm.R")
			cal_overlap(geneset_use_file)
		}
		if(bg_file){
			cat("\n*Run ORA_filter to get filtered and clustered enriched gene sets:\n\tRscript ORA_filter.R\n")
		  interest_file_id.out_LEGO_filter <- ORA_filter(interest_file_id.out_LEGO.txt, interest_file_id.out_LEGO_filter, geneset_use_file, multi)
		}else{
			cat("\n*Run ORA_filter to get filtered and clustered enriched gene sets:\n\tRscript ORA_filter.R\n")
		  interest_file_fisher_filter <- ORA_filter(interest_file_fisher, interest_file_fisher_filter, geneset_use_file, multi, 7)
		}
	}
	### finish
	cat("###### Run end:", date(), "######\n\n")
}
