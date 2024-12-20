% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gibbs_C}
\alias{gibbs_C}
\title{A Gibbs sampler using Rcpp}
\usage{
gibbs_C(N)
}
\arguments{
\item{N}{the number of samples}
}
\value{
a random matrix of nrow \code{N} and ncol 2
}
\description{
A Gibbs sampler using Rcpp
}
\examples{
\dontrun{
x<-gibbs_C(5000)
}
}
