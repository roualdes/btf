\name{tf}
\alias{tf}
\title{approximate trend filtering via MM algorithm}
\usage{
tf(y, x = NULL, k = 2, l, D = NULL, eps = 1e-08, tau = 1e-05,
  max_iter = 500)
}
\arguments{
  \item{y}{observed data}

  \item{x}{inputs corresponding to observations}

  \item{k}{order of fit}

  \item{l}{vector of penalty parameters lambda}

  \item{D}{matrix Delta of order k+1}

  \item{eps}{error adjustment to majorization function}

  \item{tau}{convergence threshold}

  \item{max_iter}{maximum number of iterations allowed}
}
\description{
Uses majorization-minimization technique to approximate
trend filtering fit.
}
\author{
Edward A. Roualdes
}

