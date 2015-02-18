\name{btf}
\alias{btf}
\title{##' Bayesian trend filtering via Eigen}
\usage{
btf(y = "vector", x = NULL, k = "int", iter = 10000,
  cond.prior = c("gdp", "dexp"), alpha = NULL, rho = NULL, D = "Matrix",
  debug = FALSE)
}
\arguments{
  \item{y}{response vector}

  \item{x}{inputs corresponding to y observations}

  \item{k}{degree of polynomial fit}

  \item{iter}{number of samples to draw from posterior}

  \item{cond.prior}{choose the conditional prior on
  f|sigma}

  \item{alpha}{shape parameter for prior on lambda}

  \item{rho}{rate parameter for prior on lambda}

  \item{D}{linear transformation of coefficients inside
  penalty}

  \item{debug}{boolean telling btf to check for NaNs or
  not}
}
\description{
Fits Bayesian trend filtering hierarchical model to
univariate function. Two conditional priors are available:
double exponential or generalized double Pareto.
}
\examples{
# Cubic trend filtering
# from genlasso::trendfilter
\dontrun{n <- 100
beta0 = numeric(100)
beta0[1:40] <- (1:40-20)^3
beta0[40:50] <- -60*(40:50-50)^2 + 60*100+20^3
beta0[50:70] <- -20*(50:70-50)^2 + 60*100+20^3
beta0[70:100] <- -1/6*(70:100-110)^3 + -1/6*40^3 + 6000
beta0 <- -beta0
beta0 <- (beta0-min(beta0))*10/diff(range(beta0))
y <- beta0 + rnorm(n)
bfit <- btf(y=y, k=3)
plot(bfit, col='grey70')}
}
\author{
Edward A. Roualdes
}
\references{
R. J. Tibshirani. Adaptive piecewise polynomial estimation
via trend filtering. The Annals of Statistics,
42(1):285-323, 2014.
}
\seealso{
\code{\link[genlasso]{trendfilter}}
}

