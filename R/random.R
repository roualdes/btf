zrnormChol <- function(N, mu, Sig) {
    ncols <- ncol(Sig)
    z <- zrnorm(N)
    mu + Sig%*%z
}

rinvGauss <- function(N, nu, lambda) {
    z <- zrnorm(N)
    z2 <- z*z
    nu2 <- nu*nu
    x <- nu + 0.5*z2*nu2/lambda - (0.5*nu/lambda)*sqrt(4.0*nu*lambda*z2+nu2*(z2*z2))
    idx <- runif(N) < nu/(nu+x)
    x[idx] <- nu2[idx]/x[idx]
    as.vector(x)
}

rinvGamma <- function(N, shape, scale) {
    1/rgamma(n = N, shape = shape, rate = scale)
}


