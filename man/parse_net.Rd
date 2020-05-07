\name{parse_net}
\alias{parse_net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert to integer ID}
\description{
  This script read in original network file and convert to integer ID.
}
\usage{
parse_net(network_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{network_file}{ a network gene interaction set object including data data.frame }
}
\value{
returns the corresponding gene interaction set Ids.
}
\author{ Zhenyu Hu }
\seealso{  \code{\link{parse_gs}} }
\examples{

    ## return a data.frame
    network <- parse_net(network_file)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
