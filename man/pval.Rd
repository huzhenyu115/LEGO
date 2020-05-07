\name{pval}
\alias{pval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Calculate a p-value}
\description{
  Calculate a p-value.
}
\usage{
pval(input, thre, met, tmp, multi)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{ a interest gene set id including data data.frame }
  \item{thre}{ a p value threshold }
  \item{met}{ a multiple testing corrections method }
  \item{tmp}{ a middle result file }
  \item{multi}{ Multi sample or not }
}
\value{
returns the corresponding p value.
}
\author{ Zhenyu Hu }
\seealso{  \code{\link{interpret_gpd}} , \code{\link{pval}} }
\examples{

    ## return a data.frame
	pval(input, thre, met, tmp, multi)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
