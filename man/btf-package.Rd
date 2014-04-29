\name{btf-package}
\alias{btf-package}
\docType{package}
\title{
  Estimate univariate non-parametric regression function with domain [0,1].
}
\description{
  Trend filtering, as developed by Ryan J. Tibshirani, estimates the function f_0 in the model: y_i = f_0(x_i) + epsilon_i, i = 1, ..., n. The btf package provides a Bayesian solution, estimating all parameters from a Gibbs sampler, thus providing confidence intervals and an estimate of the penalty parameter lambda.
}
\details{
\tabular{ll}{
Package: \tab btf\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-12-10\cr
License: \tab GPL (>= 2.0)\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
  Edward A. Roualdesp
  Maintainer: Edward A. Roualdes <edward.roualdes@uky.edu>
}
\references{
  Tibshirani, R. J. and Taylor, J. (2011), "The solution path of the
  generalized lasso", Annals of Statistics 39 (3) 1335--1371.
}
\keyword{package}
\seealso{
  \code{\link[genlasso]{trendfilter}}
}