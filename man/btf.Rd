\name{btf}
\alias{btf}
\title{Bayesian trend filtering via Eigen}
\usage{
btf(y = "vector", x = NULL, k = "int", iter = 5000,
  cond.prior = c("gdp", "dexp"), alpha = NULL, rho = NULL)
}
\arguments{
  \item{y}{response vector}

  \item{x}{inputs corresponding to y observations}

  \item{k}{degree of polynomial fit}

  \item{iter}{number of samples to draw from posterior}

  \item{cond.prior}{choose the conditional prior on
  f|sigma}

  \item{alpha}{shape parameter for prior on lambda}

  \item{rho}{rate parameter for prior on lambda}
}
\description{
need better description of this function
}
\author{
Edward A. Roualdes
}

