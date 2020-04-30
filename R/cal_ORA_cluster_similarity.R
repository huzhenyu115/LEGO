# this function aims to calculate similarity between cluster results
cal_ORA_cluster_similarity <- function(enrich_file, geneset_over_file){
#enrich_file = ARGV[0] ## output of ORA_filter.pl
#geneset_over_file = $ARGV[1]
#$geneset_over_file = "demo/GeneSet_human.txt_FC2_human_overlap_union.txt"
geneset_over_file_input <- read.table(geneset_over_file, header = F, sep = "\t", fill = TRUE)
geneset_over_file_input <- as.matrix(geneset_over_file_input)
score <- matrix(nrow=length(unique(geneset_over_file_input[,1])),ncol=length(unique(geneset_over_file_input[,2])))
row_score <- unique(geneset_over_file_input[,1])
col_score <- unique(geneset_over_file_input[,2])
for(i in length(score[,1])){
	gs1 = geneset_over_file_input[i,1]
	gs2 = geneset_over_file_input[i,2]
	score[gs1,gs2] = geneset_over_file_input[i,3]
	score[gs2,gs1] = geneset_over_file_input[i,3]
}
## enrich_file input
enrich_file_input <- read.table(enrich_file, header = F, sep = "\t", fill = TRUE)
enrich_file_input <- as.matrix(enrich_file_input)

for(i in length(enrich_file_input[,1])){
	if(str_extract(gs, "^Results for (.*)")){
		geo = str_extract(gs, "^Results for (.*)")
		next
	}
	if(str_extract(gs, "^Cluster(.*)\:")){
		cluster_id = str_extract(gs, "^Cluster(.*)\:")
		next
	}
	a = split "\t"
	for(each in a){
		tmp = split "_",each
		name = join "_",tmp[0:(length(tmp)-1)]
		result[geo,cluster_id,name] = 1
	}
}
##
all_geo = keys(result)
for(i in all_geo){
	cluster_id_1 = keys(result[[i]])
	for(j in all_geo){
		cluster_id_2 = keys(result[[j]])
		ss = 0
		for(c1 in cluster_id_1){ ## for each cluster 1, find max
			max_score = 0
			for(c2 in cluster_id_2){
				n1_arr = keys(result[i,c1])
				n2_arr = keys(result[j,c2])
				## count
				tmp <- hash()
				map{tmp{$_}=1}n1_arr
				D=0
				foreach n (n2_arr){
					if(tmp[[n]]){
						D++
					}
				}
				B = length(n1_arr)
				C = length(n2_arr)
				#A = B+C-D
				if(B<C){
					tmp_score = D/B
				}else{
					tmp_score = D/C
				}				
				t1 = 0
				for(n1 in n1_arr){
					tmp_max_score = 0
					for(n2 in n2_arr){
						if(score[n1,n2]){
							tmp_max_score = (tmp_max_score>score[n1,n2])?tmp_max_score:score[n1,n2]
						}	
						if(n1 == n2){
							tmp_max_score = 1
						}
					}
					t1 += tmp_max_score
				}
				tmp_score = t1/(length(n1_arr))
				max_score = (max_score>tmp_score)?max_score:tmp_score
			}
			#print i."\t".j."\t".max_score."\n"
			ss+=max_score
		}
		score[i,j] = ss/(length(cluster_id_1))
	}
}
## final score
for(i in all_geo){
	for(j in all_geo){
		final = (score[i,j]+score[j,i])/2
		cat(i, "\t", j, "\t", final, "\n")
	}
}
}

