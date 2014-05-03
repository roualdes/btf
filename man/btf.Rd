\name{btf}
\alias{btf}
\title{Bayesian trend filtering via Eigen}
\usage{
btf(y = "vector", x = NULL, k = "int", iter = 5000, alpha = 2,
  rho = 0.01, lambda = NULL, pb = TRUE)
}
\arguments{
  \item{y}{response vector}

  \item{x}{inputs corresponding to y observations}

  \item{k}{degree of polynomial fit}

  \item{iter}{number of samples to draw from posterior}

  \item{lambda}{set lambda to a specified values}

  \item{alpha}{shape parameter for prior on lambda}

  \item{rho}{rate parameter for prior on lambda}

  \item{pb}{logical indicating weather or not to show a
  progress bar}
}
\description{
need better description of this function
}
\author{
Edward A. Roualdes
}

