# This script read in output file, convert to original name
# usage: perl src/interpret.R <out file> <gene set ID file>
# output: 
# <out file>_txt : e.g : test_data/Int_Gene_list.txt_id.out.txt
interpret <- function(input                              ## input out file
			     	  ,id_s                              ## input gene set Id file					  
					  ,met
					  ,thre 
					  ,multi)
{
input <- input
id_s <- id_s
met <- met
thre <- thre
multi <- multi
out1 = paste(input, "_LEGO.txt", sep=".")

if(multi>0){
	int_id <- input
	int_id <- gsub("_id.out", "_id_2", int_id)         ##input.txt_id.out
	int_id_input <- read.table(file=int_id, header = F, sep = "\t", fill = TRUE)
	int_id_input <- as.matrix(int_id_input)
	toID <- hash(int_id_input[,2], int_id_input[,1])
}
## geneset_use_file_id input
id_s_input <- read.table(id_s, header = F, sep = "\t", fill = TRUE)
names(id_s_input) <- c("pathway_Name", "pathway_ID")
#id_s_input <- as.matrix(id_s_input)
#ids <- hash(id_s_input[,2], id_s_input[,1])

## convert Z to P and do adjustment
tmp = paste(input, "_tmp", sample(1:10000, 1), sep="")
cat("Rscript pval.R", input, thre, met, tmp, multi, "\n")
pval(input, thre, met, tmp, multi)
tmp_input <- read.table(file=paste(tmp, "_EdgeResults.txt", sep=""), header = T, sep = "\t", fill = TRUE)
#tmp_input <- as.matrix(tmp_input)
result <- data.frame()
for(i in 1:nrow(tmp_input)){
	if(multi > 0){
		tmp_input[i,length(tmp_input[i,])] = toID[[tmp_input[i,length(tmp_input[i,])]]]	
	}
	v <- tmp_input[i,2:length(tmp_input[i,])]
	gs <- filter(id_s_input, pathway_ID==tmp_input[i,1])[1]
	result <- rbind(result, data.frame(gs, v))
}
write.table(result, file=out1, sep = "\t", quote = FALSE, col.names = FALSE)
cat("Finish! Check output file:", out1, "\n")
system("rm -rf tmp")
EdgeResults <- paste(tmp, "EdgeResults.txt", sep="_")
system("rm -rf EdgeResults")

out <- result
return(out)

}