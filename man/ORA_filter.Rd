\name{ORA_filter}
\alias{ORA_filter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ORA filter}
\description{
  ORA filter.
}
\usage{
ORA_filter(enrich_file, output_file, geneset_file, multi = 0, p_thre = 0.1, jac_thre = 0.15, p_col = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{enrich_file}{ a enrich set object including data data.frame }
  \item{output_file}{ a output set object including data data.frame }
  \item{geneset_file}{ a gene set object including data data.frame }
  \item{multi}{ Multi sample or not }
  \item{p_thre}{ a p-value }
  \item{jac_thre}{ a jac threshold }
  \item{p_col}{ collumn number }
}
\value{
returns the filter output.
}
\author{ Zhenyu Hu }
\seealso{  \code{\link{cal_overlap_union}} }
\examples{

    ## return a data.frame
	interest_file_id.out_LEGO_filter <- ORA_filter(interest_file_id.out_LEGO.txt, interest_file_id.out_LEGO_filter, geneset_use_file, multi)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
