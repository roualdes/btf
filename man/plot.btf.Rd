\name{plot.btf}
\alias{plot.btf}
\title{plot btf obects}
\usage{
\method{plot}{btf}(x, t = NULL, burn = 1000, probs = c(0.05, 0.95),
  est = median, ...)
}
\arguments{
  \item{x}{btf object}

  \item{t}{domain of function}

  \item{burn}{size of burn-in,}

  \item{probs}{numeric 2-vector of probabilities with
  values in [0,1]}

  \item{est}{function specifying how draws from the
  posterior are summarized}

  \item{...}{extra arguments to be passed as methods}
}
\description{
plot btf obects
}
\author{
Edward A. Roualdes
}

