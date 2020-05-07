\name{pert_int}
\alias{pert_int}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ label shuffling }
\description{
  Label shuffling.
}
\usage{
pert_int(gene_id, i, perm_times)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gene_id}{ a network gene interaction set id including data data.frame }
  \item{i}{ Size }
  \item{perm_times}{ Random number }
}
\value{
returns the corresponding LEGO output.
}
\author{ Zhenyu Hu }
\seealso{  \code{\link{parse_gs}} }
\examples{

    ## return a data.frame
	tmp_input <- pert_int(gene_id, i, perm_times)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
