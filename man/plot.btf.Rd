\name{plot.btf}
\alias{plot.btf}
\title{plot btf object}
\usage{
\method{plot}{btf}(x, t = NULL, burn = 1000, est = median,
  probs = c(0.025, 0.975), ...)
}
\arguments{
  \item{x}{btf object}

  \item{t}{domain of function}

  \item{burn}{size of burn-in,}

  \item{est}{function specifying how draws from the
  posterior are summarized}

  \item{probs}{numeric 2-vector of credible interval
  probabilities; if FALSE, no credible intervals are plot}

  \item{...}{extra arguments}
}
\description{
Plots a btf object with optional credible intervals.
}
\author{
Edward A. Roualdes
}

