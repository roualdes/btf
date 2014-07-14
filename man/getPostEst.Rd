\name{getPostEst}
\alias{getPostEst}
\title{get posterior estimates from btf object}
\usage{
getPostEst(btf, parameter = c("beta", "s2", "lambda", "omega", "alpha"),
  burn = 1000, est = mean)
}
\arguments{
  \item{btf}{btf object}

  \item{parameter}{name of the paramater of interest}

  \item{burn}{number of samples to discard}

  \item{est}{estimate of the object}
}
\description{
Return posterior estimates, computed from btf posterior
draws. User can supply their own estimating function,
\code{est}, applied to the columns of \link{getPost}.
}

