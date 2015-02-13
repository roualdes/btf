\name{getPost}
\alias{getPost}
\title{get posterior draws from btf object}
\usage{
getPost(btf, parameter = c("beta", "s2", "lambda", "omega", "alpha", "rho"),
  burn = 1000)
}
\arguments{
  \item{btf}{btf object}

  \item{parameter}{name of the parameter of interest}

  \item{burn}{number of samples to discard}
}
\description{
Return posterior draws from a btf object into a coda object
}
\seealso{
\code{\link[btf]{getPostEst}}
}

