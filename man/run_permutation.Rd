\name{run_permutation}
\alias{run_permutation}
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
							, perm_times)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{interest_file_id.out}{ a network gene interaction set id including data data.frame }
  \item{geneset_use_file}{ a gene set including data data.frame }
  \item{network_file}{ a network gene interaction set including data data.frame }
  \item{geneset_file}{ a gene set including data data.frame }
  \item{main_dir}{ main directory }
  \item{perm_times}{ Random number }  

}
\value{
returns the corresponding LEGO output.
}
\author{ Zhenyu Hu }
\seealso{  \code{\link{pert_int}}, \code{\link{LEGO_noperm}}, \code{\link{summary_pert}} }
\examples{

    ## return a data.frame
	mid_prefix = "mid_data/"
	interest_file_id.out <- run_permutation(interest_file_id.out, geneset_use_file, network_file, geneset_file, main_dir, perm_times)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
