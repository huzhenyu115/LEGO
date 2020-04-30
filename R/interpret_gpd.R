#!/usr/bin/perl -w
# This script read in output file, convert to original name
# usage: perl src/interpret.pl <out file> <gene set ID file>
# output: 
# <out file>_txt : e.g : test_data/Int_Gene_list.txt_id.out.txt
interpret_gpd <- function(input_id                              ## input out file
			     	      ,each_gs                              ## input gene set Id file
					      ,met
					      ,thre 
					      ,multi
						  ,mid_prefix)
{
input_id <- input_id
each_gs <- each_gs
id_s <- paste(each_gs, "_id", sep="") 
gene_id <- paste(each_gs, "_gene_id", sep="")
input <- paste(input_id, ".out", sep="")
out1 <- paste(input, "_LEGO.txt", sep="")
met <- met
thre <- thre
multi <- multi
mid_prefix <- mid_prefix

if(multi>0){
	int_id_input <- interest_gene$interest_gene_sample
	int_id_input <- as.matrix(int_id_input)
	toID <- hash(int_id_input[,2], int_id_input[,1])
}
## id_s input
id_s_input <- read.table(file=id_s, header = F, sep = "\t", fill = TRUE)
#names(id_s_input) <- c("pathway_Name", "pathway_ID")
#id_s_input <- as.matrix(id_s_input)
ids <- hash(id_s_input[,2], id_s_input[,1])
gs_id_num <- length(ids)
## convert Z to P and do adjustment
GeneSet_des_input <- read.table("GeneSet_des.txt", header = F, sep = "\t", fill = TRUE)
GeneSet_des_input <- as.matrix(GeneSet_des_input)
if(length(GeneSet_des_input) < 3){next}
des <- hash(GeneSet_des_input[,1], paste(GeneSet_des_input[,2], GeneSet_des_input[,3], sep="_"))
names(GeneSet_des_input) <- c("GO_ID", "GO_Name", "Class")
##
tmp <- paste(input, "_tmp", sample(1:10000, 1), sep="")
cat("Rscript pval_gpd.R", input, thre, met, tmp, multi, mid_prefix, "\n")
pval_gpd(input, thre, met, tmp, multi, mid_prefix)
########################## transfer id to gene
use_input1 <- read.table(file=input_id, header = F, sep = "\t", fill = TRUE)     ## input gene id
names(GeneSet_des_input) <- c("gene_ID")
use_input <- paste(use_input1, collapse=";")
#use_input_hash <- hash(use_input1[,1], rep(1, times=length(use_input1[,1])))    ## use input gene id 

use_input_num <- nrow(use_input1) 
## gene_id input
gene_id_input <- read.table(file=gene_id, header = F, sep = "\t", fill = TRUE)     
#names(id2gene) <- c("gene_Name", "gene_ID")
id2gene <- hash(gene_id_input[,2], gene_id_input[,1])    
gene_id_num <- length(id2gene) 
## ach_gs_id_gs input
each_gs_id_gs_input <- paste(each_gs, "id_gs", sep="_")
belong <- read.table(file=each_gs_id_gs_input, header = F, sep = "\t", fill = TRUE)
names(belong) <- c("gene_ID", "pathway_ID")

######################### get results
tmp_out <- paste(tmp, "EdgeResults.txt", sep="_")
tmp_out_input <- read.table(file=tmp_out, header = T, sep = "\t", fill = TRUE)
ori_gs_arr <- tmp_out_input$GeneSet_des
use_gs_num <- nrow(tmp_out_input)
use_gs <- paste(ori_gs_arr, collapse=";")
##
system("extract_CS paste(each_gs, "id_gs_CS", sep="_") gs_id_num use_gs_num use_gs gene_id_num use_input_num use_input") ## neighbor: gsid\tgeneid\tweight
txt <- read.table(file=paste(each_gs, "id_gs_CS", sep="_"), header = F, sep = "\t", fill = TRUE)
names(txt) <- c("gene_ID", "pathway_ID", "weight")

##output
out <- data.frame()
for(g in 1:nrow(tmp_out_input)){
	x <- filter(txt, pathway_ID==tmp_out_input[g,1])[1]
	y <- filter(belong, pathway_ID==tmp_out_input[g,1])[1]
	same_gene <- intersect(x, y)
	diff_gene <- setdiff(x, y)
	## score_ov
	score_ov <- c()
	score_nb <- c()
	for(i in same_gene){
	score_ov <- paste(score_ov, id2gene[[i]], collapse=",") 
	}
	## score_nb
	for(j in diff_gene){
	score_nb <- paste(score_nb, id2gene[[j]], collapse=",") 
	}
	score <- cbind(score_ov, score_nb)
	if(multi > 0){
		tmp_out_input[g,ncol(tmp_out_input)] = toID[[tmp_out_input[g,ncol(tmp_out_input)]]]	
	}
	v = tmp_out_input[g, 2:ncol(tmp_out_input)]
	ori_gs = tmp_out_input[g,1]
	gs = ids[[ori_gs]]
	if(des[[gs]]){gs=paste(gs, des[[gs]], sep="_")}
	out <- rbind(out, data.frame(gs, v, score_ov, score_nb)) 
 
}
if(multi > 0){
	colnames(out) <- c("GeneSet_des", "Ori_p", "Adj_p",	"gsSize", "gsSize_net",	"intSize", "intSize_net", "OverlapSize", "InterestingGeneListID", "Essential overlapped genes (sorted by weight)", "Essential neighbor genes (sorted by weight)")
}else{
	colnames(out) <- c("GeneSet_des", "Ori_p", "Adj_p",	"gsSize", "gsSize_net",	"intSize", "intSize_net", "OverlapSize", "Essential overlapped genes (sorted by weight)", "Essential neighbor genes (sorted by weight)")
}

write.table(out, file=out1, sep = "\t", quote = FALSE, col.names = FALSE)
##
cat("Finish! Check output file:", out1, "\n")
system("rm -rf tmp")
EdgeResults <- paste(tmp, "EdgeResults.txt", sep="_")
system("rm -rf EdgeResults")

return(out)

}