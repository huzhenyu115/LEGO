\name{parse_gs}
\alias{parse_gs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert to integer ID}
\description{
  This script read in original gene set file and convert to integer ID.
}
\usage{
parse_gs(geneset_file, network_file_id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneset_file}{ a network gene set object including data data.frame }
  \item{network_file_id}{ a network gene interaction set id object including data data.frame }
}
\value{
returns the corresponding gene set Ids.
}
\author{ Zhenyu Hu }
\seealso{  \code{\link{parse_net}} }
\examples{

    ## return a data.frame
    geneset <- parse_gs(geneset_file, network$network_file_id)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
