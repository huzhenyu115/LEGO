# LEGO
Gene set function enrichment
* INTRODUCTION
----------------------------------------------------------------------

We have developed a novel method named LEGO (functional Link Enrichment of Gene Ontology or gene sets). Incorporating a network-based gene-weighting scheme, LEGO measures the overlaps between the interesting genes and a given gene set by integrating both the gene weights and gene-gene interactions.   

Please Cite:

LEGO: a novel method for gene set over-representation analysis by incorporating network-based gene weights. 

* HOWTO COMPILE THE SOURCE AND INSTALL	
----------------------------------------------------------------------

	You will need a Unix environment with a working C++ compiler to
	compile the source files in src/. We recommend GCC>=4.2.
	Besides, Perl and R are necessary for running this program.

* USAGE
----------------------------------------------------------------------

In the main directory, you could type:
  LEGO(<network file> <geneset file> <interest file> [options])
and you will see the help message.

The usage for this program: 

###############################################################################################
	
  LEGO(<network file> <geneset file> <interest file> [options])

	Options:  

	pre: <whether or not to pre-run to generate mid files,1-yes,0-no(e.g:0)> 
	multi: <multi or not:1-yes,0-no(e.g:0)> 
	bgNE: <bgNE for the network,e.g:0.25> 
	min: <minimum gene set size,e.g:5> 
	max: <maximum gene set size,e.g:10000> 
	adj: <adjusted methods,e.g:fdr> 
	p_thre: <p value cutoff,e.g:0.05> 
	fisher: <whether or not run fisher exact test, 1-yes,0-no,e.g: 0> 
	filter: <whether or not run result cluster and filter step, 1-yes, 0-no, e.g: 1> 
	bg_file: <background file, e.g: bg.txt, if do not have background file, leave it blank>

###############################################################################################

And we also provide a version that do not use permutation strategy:

The usage for this program:

###############################################################################################

LEGO_noperm(<network file> <geneset file> <interest file> [options])

        Options:

        pre: <whether or not to pre-run to generate mid files,1-yes,0-no(e.g:0)>
        multi: <multi or not:1-yes,0-no(e.g:0)>
        bgNE: <bgNE for the network,e.g:0.25>
        min: <minimum gene set size,e.g:5>
        max: <maximum gene set size,e.g:10000>
        adj: <adjusted methods,e.g:fdr>
        p_thre: <p value cutoff,e.g:0.05>
        fisher: <whether or not run fisher exact test, 1-yes,0-no,e.g: 0>
        filter: <whether or not run result cluster and filter step, 1-yes, 0-no, e.g: 0>
        bg_file: <background file, e.g: bg.txt, if do not have background file, leave it blank> \n";

###############################################################################################

	
Example1 (for single interesting gene list file): 
LEGO("data/FC2_human.txt", "data/GeneSet_human.txt", "data/input.txt")

Example2 (for multiple interesting gene list file): 
LEGO("data/FC2_human.txt", "data/GeneSet_human.txt", "data/input_mul.txt", multi=1)

Example3 (for single interesting gene list file and filtered results): 
LEGO("data/FC2_human.txt", "data/GeneSet_human.txt", "data/input.txt", filter=1)

Example4 (for multiple interesting gene list file and filtered results): 
LEGO("data/FC2_human.txt", "data/GeneSet_human.txt", "data/input_mul.txt", multi=1, filter=1)

Example5 (for single interesting gene list file + bg file): 
LEGO("data/FC2_human.txt", "data/GeneSet_human.txt", "data/input.txt", bg_file="demo/bg.txt")

Example6 (for single interesting gene list file with background list): 
LEGO_noperm("data/FC2_human.txt", "data/GeneSet_human.txt", "data/input.txt", bg_file="demo/bg.txt)

Example7 (for multiple interesting gene list file with background list): 
LEGO_noperm("data/FC2_human.txt", "data/GeneSet_human.txt", "data/input_mul.txt", multi=1, bg_file="data/bg.txt")

Tips: 
	1. Output files: <interest file>_LEGO.txt, if you choose to filter & cluster enriched gene sets, another two files will be generated: <interest file>_LEGO_filter.txt <interest file>_LEGO_filter_cluster.txt (set -filter 1)
	2. The first time to run this program will take about several minutes to generate mid-files(according to the size of the network and genesets). If you want to re-generate the mid files, set `-pre_run to 1`; 
	If you have already generated the mid-files for network and genesets, it will only take several seconds to calculate a new input gene list (set `-pre_run to generate mid files` to 0)
	3. The first time to use LEGO(), it will take maybe several minutes to generate the permutation results. We provide the pre-generated mid files for data/GeneSet_human.txt. 
	4. If you want to use the background file, you could only use LEGO_noperm();
	5. The gene ID in the network, gene set and input gene list must be the same ID system.
