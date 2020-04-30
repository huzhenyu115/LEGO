# Usage:  enrich.pl  <database file> <input_file_name> <output_file_name>  [options]
# database file :  database file for enrichment analyais
# input_file_name:   input file name
# output_file_name: 	output file name
# Options:
# item   (item column number,default = 1)
# des    (des column number, default = 2)
# des_file (files to save description for genesets)
# bg     (bacground_file_name,default=no background file)
# list   (1-y|2-n,default=n)
# thre   (threshold,default = 0.01)
# max    (maximun size for des, default = no limitation)
# min    (minmun size for des, default = no limitation)
# adj    (need adj p value or not ,0-no,1-BH,2-FDR,default = no)
# low    (low existence for match,default = 0)
# multi  (whether the input has multi groups, default=0)
# filter ( use what filter or not, default = 0 (do not use))
# filter_column (filter column)
enrich <- function(database_file
				           ,input_file
                   ,item_col    = 1
                   ,des_col     = 2
                   ,des_file    = ""
                   ,list        = 2
                   ,thre        = 0.01
                   ,background  = "no"
                   ,min_size    = 0
                   ,max_size    = 10000000
                   ,need_adj    = "none"
                   ,low_match   = 0
                   ,filter      = 0
                   ,filter_col  = 0
                   ,multi   = 0)
{
	library(stringr)
	library(hash)
	#### method type
	method_all <- c("fdr","holm","hochberg","hommel","BH","BY")
	use_method_all <- hash(method_all, rep(1,times=length(method_all)))
  if(has.key(need_adj, use_method_all)){
		adj_method <- need_adj
	}else{
		adj_method <- "none"
	}
	#adj_method

	####### read in background (if $background != 0)
	if(background != "no"){
		background_input <- read.table(file=background, header = F, sep = "\t", fill = TRUE)
		background_input <- as.matrix(background_input)
		exist_back <- hash(toupper(background_input[,1]), rep(1, times=length(background_input[,1])))
		cat("Background items:", length(exist_back) ,"\n")
	}else{
		cat("no background\n")
	}
	###### read in des file
	if(des_file != ""){
		des_file_input <- read.table(file=des_file, header = F, sep = "\t", fill = TRUE)
		des_file_input <- as.matrix(des_file_input)
		des <- hash()
		for(i in 1:length(des_file_input)){
			if(length(des_file_input[i,]) < 3){
				next
			}else{
				.set(des, keys = a[1], values = paste(a[2], a[3], sep="_"))
			}
		}
	}
	###### a: whole gene numbers
	#database_file_input <- read.table(database_file, header = F, sep = "\t", fill = TRUE)
	database_file_input <- database_file
	database_file_input <- as.matrix(database_file_input)
	exist <- hash()
	all <- hash()
	all_size <- hash()

	for(i in 1:length(database_file_input[,1])){
		item <- toupper(database_file_input[i,item_col])
		des_p <- database_file_input[i,des_col]
		if(item == "."){next}
		if(des_file != ""){
			if(has.key(des_p, des)){des = paste(des_p, des[[des_p]], sep="_")}
		}
		if(filter){
			use_filter <- toupper(database_file_input[filter_col+length(database_file_input[i,])])
			if(filter != use_filter){
				next
			}
		}
		if(grepl(" ", des_p)){next}
		if(background != "no"){
			if(hash.key(item, exist_back)){
				.set(exist, keys = item, values = 1)
				if(!grepl("\\;", des_p)){
					if(has.key(des_p, all)){
						all[[des_p]] <- c(all[[des_p]], item)
						.set(all, keys = des_p, values = all[[des_p]])
					}else{
						.set(all, keys = des_p, values = item)
					}
				}else{
					des_a <- strsplit(des_p, ";")
					for(m in des_a){
						if(has.key(m, all)){
							all[[m]] <- c(all[[m]], item)
							.set(all, keys = m, values = all[[m]])
						}else{
							.set(all, keys = m, values = item)
						}
					}
				}
			}
			if(!grepl("\\;", des_p)){
				if(has.key(des_p, all_size)){
					all_size[[des_p]] <- c(all_size[[des_p]], item)
					.set(all_size, keys = des_p, values = all_size[[des_p]])
				}else{
					.set(all_size, keys = des_p, values = item)
					}
			}else{
				des_a <- strsplit(des_p, ";")
				for(n in des_a){
					if(has.key(n, all_size)){
						all_size[[n]] <- c(all_size[[n]], item)
						.set(all_size, keys = n, values = all_size[[n]])
					}else{
						.set(all_size, keys = n, values = item)
					}
				}
			}
		}
		if(background == "no"){
			.set(exist, keys = item, values = 1)
			if(!grepl("\\;", des_p)){
				if(has.key(des_p, all)){
					all[[des_p]] <- c(all[[des_p]], item)
					.set(all, keys = des_p, values = all[[des_p]])
				}else{
					.set(all, keys = des_p, values = item)
				}
			}else{
				des_a <- strsplit(des_p, ";")
				for(k in des_a){
					if(has.key(k, all)){
						all[[k]] <- c(all[[k]], item)
						.set(all, keys = k, values = all[[k]])
					}else{
						.set(all, keys = k, values = item)
					}

				}
			}
		}
	}
	a <- length(keys(exist))
	cat("All items for analysis:", a, "\n")
	if(background == "no"){
		all_size <- all
	}

	no <- hash()
	for(each_des in keys(all_size)){
		b <- length(all_size[[each_des]])
		if(b < min_size || b > max_size){
			.set(no, keys = each_des, values = 1)
		}
	}
	#######input list c--> whole gene numbers in the list
	#insterest_input <- read.table(input_file, header = F, sep = "\t", fill = TRUE)
	#insterest_input <- as.matrix(insterest_input)
	insterest_input <- interest_file
	all_count <- data.frame()
	for(i in 1:nrow(insterest_input)){
		gene <- toupper(insterest_input[i,1])
		if(multi){
			mark <- unique(insterest_input[,2])
			group <- insterest_input[i,2]
			if(background != "no"){
				if(has.key(gene, exist_back)){
					all_count <- rbind(all_count, data.frame(group=group, gene=gene))
				}
			}else{
				all_count <- rbind(all_count, data.frame(group=group, gene=gene))

			}
		}else{
			mark <- 1
			if(background != "no"){
				if(has.key(gene, exist_back)){
					all_count <- rbind(all_count, data.frame(group=1, gene=gene))
				}
			}else{
					all_count <- rbind(all_count, data.frame(group=1, gene=gene))

			}
		}
	}
	##
	ng <- mark
	cat("Number of groups:", ng, "\n")

	##
	A = a
	if(background != "no"){
		A = length(exist_back)
		if(!multi){C = mark}
	}


	## for each group
	result <- data.frame()
	for(group in unique(all_count$group)){      ## for each group
		j <- 0
		for(each_des in keys(all)){
			if(has.key(each_des, no)){next}
			b <- length(all[[each_des]])
			gene_list <- filter(all_count, group==group)[2]
			c <- nrow(gene_list)
			if(background != "no"){
				if(multi){
					c <- length(gene_list)
				}else{
					c <- C
				}
			}
			gene_all <- intersect(as.matrix(as.data.frame(all[[each_des]])), as.matrix(gene_list))
			d <- length(gene_all)
			if(d < low_match){
				next
			}
			if(d != 0){
				gene_all_list <- paste(gene_all, collapse=",")
			}else{
				gene_all_list <- NA
			}

			##
			a <- A
			input <- c(a,b,c,d)
			input <- as.matrix(input)
			p <- apply(input,2,function(x){
					fisher.test(matrix(x, nrow=2),alternative="greater")$p.value
				})
			adj.p <- p.adjust(p,method=adj_method)
			ratio <- (d / c)/(b / a)
			result <- rbind(result, data.frame(each_des, a, b, c, d, p, adj.p, ratio, gene_all_list))
			j <- j+1
		}

		if(j == 0){
			cat("No enriched ! Change the parameters & try!\n")
			next
		}
	}

	#########enrichment --> b-> gene numbers in a GO / d-> gene numbers in a GO in the list
	if(multi){
		colnames(result) <- c("Name", "Total_Item", "Num_item", "Num_list", "Num_list_item", "Ori_p", "Adj_p", "Enrich_score", "Gene_list", "Group_ID")
	}else{
		colnames(result) <- c("Name", "Total_Item", "Num_item", "Num_list", "Num_list_item", "Ori_p", "Adj_p", "Enrich_score", "Gene_list")
	}

	#output
	out <- result
	return(out)
}
