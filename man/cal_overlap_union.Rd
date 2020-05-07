\name{cal_overlap_union}
\alias{cal_overlap_union}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate correlation between pathways }
\description{
  This script calculate correlation between pathway and pathway.
}
\usage{
cal_overlap_union(geneset_use_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneset_use_file}{ a gene set object including data data.frame }
}
\value{
returns the correlation coefficient between pathways.
}
\author{ Zhenyu Hu }
\seealso{  \code{\link{parse_gs}} }
\examples{

    ## return a data.frame
    if(!exists("geneset_use_file_overlap_union")){
		geneset_use_file_overlap_union <- cal_overlap_union(geneset_use_file)
	}
	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
