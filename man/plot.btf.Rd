\name{plot.btf}
\alias{plot.btf}
\title{plot btf obects}
\usage{
\method{plot}{btf}(x, y, chain = 1, burn = 100, probs = c(0.025, 0.975),
  est = mean)
}
\arguments{
  \item{x}{coda object that contains a btf history}

  \item{y}{response vector for btf model}

  \item{chain}{number of chains included in x}

  \item{burn}{size of burn-in,}

  \item{probs}{numeric 2-vector of probabilities with
  values in [0,1]}
}
\description{
plot btf obects
}
\author{
Edward A. Roualdes
}

