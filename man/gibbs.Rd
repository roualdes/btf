\name{gibbs}
\alias{gibbs}
\title{Gibbs sampler framework}
\usage{
gibbs(M, fn, individual, toCODA)
}
\arguments{
  \item{M}{number of iterations}

  \item{fn}{function to update parameters in-place}

  \item{individual}{member of the population who's state is
  updated by fn}

  \item{toCODA}{a function to convert state history to CODA
  object}
}
\description{
Gibbs sampler framework
}
\author{
Edward A. Roualdes
}

