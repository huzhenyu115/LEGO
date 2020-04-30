######GeneSet_des.txt file input
add_des <- function(){
GeneSet_des <- read.table(file="GeneSet_des.txt", header = F, sep = "\t", fill = TRUE)
GeneSet_des <- as.matrix(GeneSet_des)
des <- hash()
for(i in length(GeneSet_des[,1])){
	if(length(GeneSet_des[i,1]) < 2){next}
	.set(des, keys = GeneSet_des[i,1], values = paste(GeneSet_des[i,2], GeneSet_des[i,3], sep="_"))
}

a <- read.table(file=ARGV[0], header = F, sep = "\t", fill = TRUE)
a <- as.matrix(a)
for(i in 1:length(a)){
	#$pv = $a[3];
	#if($pv>0.05){next;}
	gs = a[i,1]
	# Cluster419_(265):GO:0042634_0.095617_
	if(str_extract(gs, "^GO:.......")){ 
		gs = paste(gs, des[[gs]], sep="_")
		a[i,1] = gs
	}
	n = paste(a[i,], sep="\t")
	cat(n, "\n")
}
}
