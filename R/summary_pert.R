# Rscript R/summary_pert.R each tmp_input final_file
summary_pert <- function(tmp_input_id.out, final_file, main_dir){
	##
	library(multicore);
	file_perm <- tmp_input_id.out;
	out_file <- final_file;
	path <- main_dir;
	input_info <- as.numeric(./src/extract_info(file_perm));
	all_gs <- as.numeric(./stc/extract_info_gs(file_perm));
	save(input_info,all_gs,file=out_file);
	## output
}
