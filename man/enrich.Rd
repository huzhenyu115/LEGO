\name{enrich}
\alias{enrich}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gene set enrichment analysis}
\description{
  This script performs gene set enrichment analysis.
}
\usage{
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
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{database_file}{ a gene set object including data data.frame }
  \item{input_file}{ a interest gene  set object including data data.frame }
  
}
\value{
returns the corresponding gene set enrichment output.
}
\author{ Zhenyu Hu }

\examples{

    ## return a data.frame
	if(fisher){
		interest_file_fisher <- enrich(database_file=geneset_use_file, input_file=interest_file)
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
