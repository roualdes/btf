\name{gibbs}
\alias{gibbs}
\title{Gibbs sampler framework}
\usage{
gibbs(M, fn, state, toCODA, pb = TRUE)
}
\arguments{
  \item{M}{number of iterations}

  \item{fn}{function to update parameters in-place}

  \item{state}{member of the population who's state is
  updated by fn}

  \item{toCODA}{a function to convert state history to CODA
  object}

  \item{pb}{logical indicating weather or not to show a
  progress bar}
}
\description{
implements Gibbs iterations based on updating function fn
}
\author{
Edward A. Roualdes
}

