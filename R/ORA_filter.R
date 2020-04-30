## default
ORA_filter <- function(enrich_file
					   ,output_file
					   ,geneset_file
					   ,multi = 0
					   ,p_thre = 0.1
					   ,jac_thre = 0.15
					   ,p_col = 2)
{
###########
geneset_GSM_overlap <- paste("geneset_file", "overlap_union.txt", sep="_") ## BIOCARTA_RELA_PATHWAY   BIOCARTA_ATM_PATHWAY    Gene_set_overlap    0.18
geneset_module_file <- paste(geneset_GSM_overlap, jac_thre, "result_iNP_MSG", sep="_")
if(!file.exists(geneset_GSM_overlap)){
	cat("Rscript cal_overlap_union.r", geneset_file, "\n")
	cal_overlap_union(geneset_file)
}
if(!file.exists(geneset_module_file)){
	tmp = paste("tmp", sample(1:10000, 1), sep="")
	out_tmp = paste(tmp, "result_iNP_MSG", sep="_")
	GSM_input <- read.table(file=geneset_GSM_overlap, header = F, sep = "\t", fill = TRUE)
	GSM_input <- as.matrix(GSM_input)
	out <- GSM_input[which(GSM_input[,3] >=jac_thre),]
	write.table(out, file=tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
	cat("iNP", tmp, "1 0 tmp\n")
	system("iNP tmp 1 0 tmp")
	system("mv out_tmp geneset_module_file")
	system("rm -rf tmp*")
}
## read in module file
all_gs <- read.table(file=geneset_file, header = F, sep = "\t", fill = TRUE)
all_gs <- as.matrix(all_gs)
m <- 1
geneset_module <- read.table(file=geneset_module_file, header = F, sep = "\t", fill = TRUE)
geneset_module <- as.matrix(geneset_module)
cat("Modularity:",geneset_module[length(geneset_module[,1]),1])

belong <- hash()
m <- 1
for(i in 1:nrow(geneset_module)){
	for(j in 1:ncol(geneset_module)){
		if(geneset_module[i,j] != ""){
			.set(belong, keys=geneset_module[i,j], values=m)
			m <- m+1
		}
	}
}

n1 <- m
cat("There are", n1, "gene sets in iNP results\n")
for(gs in all_gs[,2]){
	if(has.key(gs, belong)){
		next
	}
	m <- m+1
	.set(belong, keys = gs, values = m)
}

##
cat("In total", m, "modules\n")
## install.packages("stringr")
library(stringr)
enrich_input <- read.table(file=enrich_file, header = T, sep = "\t", fill = TRUE)
enrich_input <- as.matrix(enrich_input)
all_int <- c()
tmp1 <- hash()
final <- list()
count <- 1
for(i in 1:nrow(enrich_input)){
	if(multi == 1){
		int <- enrich_input[i,nrow(enrich_input)]
	}else{
		int <- "ONE"
	}
	if(!has.key(int, tmp1)){
		all_int <- c(all_int, int)
		.set(tmp1, keys = int, values = 1)
	}
	ori_gs <- enrich_input[i,1]
	gs <- enrich_input[i,1]
	if(grepl("GO:.......", ori_gs)){
		gs <- grepl("GO:.......", ori_gs)
	}
	if(p_col == 0){
		pv <- 1
	}else{
		pv <- as.numeric(enrich_input[i,p_col])
	}
	if(abs(pv) <= p_thre){
		bgs <- belong[[gs]]
		if(multi == 1){
			final[[count]] <- int
			final[[count]][2] <- bgs
			final[[count]][3] <- ori_gs
			final[[count]][4] <- pv
		}else{
			final[[count]] <- int
			final[[count]][2] <- bgs
			final[[count]][3] <- ori_gs
			final[[count]][4] <- pv
		}
		count <- count+1
	}
}
##
mark <- hash()
mark2 <- hash()
for(i in 1:length(final)){
	if(has.key(final[[i]][2],mark)){
		a <- strsplit(mark[[final[[i]][2]]], ":")[[1]][2]
		if(a > final[[i]][4]){
			.set(mark, keys=final[[i]][2], values=paste(final[[i]][3],final[[i]][4], sep=":"))
			.set(mark2, keys=final[[i]][2], values=paste(mark[[final[[i]][2]]], paste(final[[i]][3],final[[i]][4], sep="_"), sep=";"))
		}
		.set(mark2, keys=final[[i]][2], values=paste(mark[[final[[i]][2]]], paste(final[[i]][3],final[[i]][4], sep="_"), sep=";"))
	}else{
		.set(mark, keys=final[[i]][2], values=paste(final[[i]][3],final[[i]][4], sep=":"))
		.set(mark2, keys=final[[i]][2], values=paste(final[[i]][3],final[[i]][4], sep="_"))
	}
}
mark1 <- hash()
for(i in keys(mark)){
	.set(mark1, keys=i, values=strsplit(mark[[i]], ":")[[1]][1])
}

output_out <- data.frame()
for(int in all_int){
	for(i in 1:length(final)){
		mark_gs <- mark1[[final[[i]][2]]]
		n <- ""
		bgs <- final[[i]][2]
		gs <- final[[i]][3]
		pv <- final[[i]][4]
		n <- paste(n, paste(gs, pv, sep="_"), sep=";")
		if(multi == 1){
				output_out = rbind(output_out, data.frame(gs, pv, bgs, mark_gs, int))
		}else{
				output_out = rbind(output_out, data.frame(gs, pv, bgs, mark_gs))
		}
		tmp = paste(mark_gs, "_", strsplit(mark[[final[[i]][2]]],":")[[1]][2], "\n", n, sep="")
	}
}
k <- 0
cluster_out <- list()
for(int in all_int){
	for(i in keys(mark2)){
		k <- k+1
		if(multi == 1){
			#cluster_out = c(tab_out, paste("Result for", paste(int, ":", sep=""), sep=" "))
		}
		cluster_out[[k]] = paste(paste("Cluster", k, sep=""),mark[[i]],sep=":" )
		cluster_out[[k]][2] = mark2[[i]]
	}
}
write.table(output_out, file=output_file, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(cluster_out, file=paste(output_file, "tab.txt", sep="_"), sep = "\t", quote = FALSE, col.names = FALSE)

out <- list(Filter=output_out, Filter_tab=cluster_out)
return(out)

}
