\name{getPostEst}
\alias{getPostEst}
\title{get posterior estimates from btf object}
\usage{
getPostEst(btf, parameter = c("beta", "s2", "lambda", "omega", "alpha",
  "rho"), burn = 1000, est = NULL)
}
\arguments{
  \item{btf}{btf object}

  \item{parameter}{name of the paramater of interest}

  \item{burn}{number of samples to discard}

  \item{est}{function to estimate from the posterior
  samples of \code{parameter}}
}
\description{
Return posterior estimates, computed from btf posterior
draws. User can supply their own estimating function,
\code{est}, applied to the columns of
\code{\link[btf]{getPost}}.
}
\seealso{
\code{\link[btf]{getPost}}
}

