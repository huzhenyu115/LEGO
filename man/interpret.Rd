\name{interpret}
\alias{interpret}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert to original name}
\description{
  This script read in output file, convert to original name.
}
\usage{
interpret(input, id_s, met, thre, multi)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{ a interest gene set id including data data.frame }
  \item{id_s}{ a gene set object including data data.frame }
  \item{met}{ a multiple testing corrections method }
  \item{thre}{ a p value threshold }
  \item{multi}{ Multi sample or not }
}
\value{
returns the corresponding original gene set name.
}
\author{ Zhenyu Hu }
\seealso{  \code{\link{interpret_gpd}} , \code{\link{pval}} }
\examples{

    ## return a data.frame
    input_LEGO <- interpret(interest_file_id,  geneset_use_file, adj_meth, p_thre, multi)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
