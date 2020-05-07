\name{parse_int_multi}
\alias{parse_int_multi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert to integer ID}
\description{
  This script read in interesting gene list file and convert to integer ID.
}
\usage{
parse_int_multi(interest, geneset_gene_id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{interest}{ a interest gene set object with multi-sample including data data.frame }
  \item{geneset_gene_id}{ a interest geneList Id object including data data.frame }
}
\value{
returns the corresponding gene set Ids.
}
\author{ Zhenyu Hu }
\seealso{  \code{\link{parse_int_multi_bg}} }
\examples{

    ## return a data.frame
    interest_gene_id <- parse_int_multi(interest_file, geneset$geneset_use_file_gene_id)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
