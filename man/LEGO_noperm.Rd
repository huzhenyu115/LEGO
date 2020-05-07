\name{LEGO_noperm}
\alias{LEGO_noperm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ noperm LEGO}
\description{
  This script run noperm LEGO.
}
\usage{
LEGO_noperm(network_file
				 ,geneset_file
                 ,interest_file
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
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{network_file}{ a network gene interaction set object including data data.frame }
  \item{geneset_file}{ a gene set object including data data.frame }
  \item{interest_file}{ a interest gene set object including data data.frame }
  \item{pre_run}{ whether or not to pre-run to generate mid files,1-yes,0-no(e.g:0) }
  \item{min}{ minimum gene set size,e.g:5 }
  \item{max}{ maximum gene set size,e.g:10000 }
  \item{multi}{ multi or not:1-yes,0-no(e.g:0) }
  \item{bgNE}{ bgNE for the network,e.g:0.25 }
  \item{adj_meth}{ adjusted methods,e.g:fdr }
  \item{p_thre}{ p value cutoff,e.g:0.05 }
  \item{fisher}{ whether or not run fisher exact test, 1-yes,0-no,e.g: 0 }
  \item{filter}{ whether or not run result cluster and filter step, 1-yes, 0-no, e.g: 0 }
  \item{bg_file}{ background file, e.g: bg.txt, if do not have background file, leave it blank }
  \item{perm_times}{ random number }
  \item{no_R}{ label }
  
}
\value{
returns the corresponding LEGO_noperm output.
}
\author{ Zhenyu Hu }
\seealso{  \code{\link{LEGO}} }
\examples{

    ## return a data.frame
	LEGO_noperm(net_file, geneset_file, tmp_input, multi=1, noR=0, min=0, max=100000)
	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
