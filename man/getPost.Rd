\name{getPost}
\alias{getPost}
\title{extract posterior samples from btf object}
\usage{
getPost(btf, parameter = c("beta", "s2", "lambda2", "omega2"), burn = 1000)
}
\arguments{
  \item{btf}{object}

  \item{parameter}{name of the parameter of interest}

  \item{burn}{number of samples to discard}
}
\description{
extract posterior samples from btf object
}

