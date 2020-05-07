\name{summary_pert}
\alias{summary_pert}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract information about gene interaction weight and pathway ids}
\description{
  Extract information about gene interaction weight and pathway ids.
}
\usage{
summary_pert(tmp_input_id.out, final_file, main_dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tmp_input_id.out}{ a middle file including data data.frame }
  \item{final_file}{ a result output  }
  \item{main_dir}{ a main work directory }

}
\value{
returns the corresponding gene interaction weight and pathway ids.
}
\author{ Zhenyu Hu }

\examples{

    ## return a data.frame
    summary_pert(tmp_input_id.out, final_file, main_dir)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
