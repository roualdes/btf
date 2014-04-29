\name{btf}
\alias{btf}
\title{Bayesian trend filtering via Eigen}
\usage{
btf(y = "vector", x = NULL, k = "int", iters = 5000, rho = 1,
  delta = 0, lambda = NULL, pb = TRUE)
}
\arguments{
  \item{y}{response vector}

  \item{x}{inputs corresponding to y observations}

  \item{k}{degree of polynomial fit}

  \item{iters}{number of samples to draw from posterior}

  \item{lambda}{set lambda to a specified values}

  \item{rho}{shape parameter for prior on lambda}

  \item{delta}{rate parameter for prior on lambda}

  \item{pb}{logical indicating weather or not to show a
  progress bar}
}
\description{
need better description of this function
}
\author{
Edward A. Roualdes
}

