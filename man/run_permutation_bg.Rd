\name{run_permutation_bg}
\alias{run_permutation_bg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ label shuffling }
\description{
  Label shuffling.
}
\usage{
run_permutation <- function(interest_file_id.out
							, geneset_use_file
							, network_file
							, geneset_file
							, main_dir
							, bg_file
							, perm_times)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{interest_file_id.out}{ a network gene interaction set id including data data.frame }
  \item{geneset_use_file}{ a gene set including data data.frame }
  \item{network_file}{ a network gene interaction set including data data.frame }
  \item{geneset_file}{ a gene set including data data.frame }
  \item{main_dir}{ main directory }
  \item{bg_file}{ a background file }  
  \item{perm_times}{ Random number }  

}
\value{
returns the corresponding LEGO output.
}
\author{ Zhenyu Hu }
\seealso{  \code{\link{pert_int}}, \code{\link{LEGO_noperm}}, \code{\link{summary_pert}} }
\examples{

    ## return a data.frame
	bg = strsplit(bg_file, "/")
	mid_prefix = paste("mid_data/", bg[[1]][2], "_", sep= "_")
	interest_file_id.out <- run_permutation_bg(interest_file_id.out, geneset_use_file, network_file, geneset_file, main_dir, bg_file, perm_times)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
