##' Bayesian trend filtering via Eigen
##'
##' need better description of this function
##'
##' @param y response vector
##' @param x inputs corresponding to y observations
##' @param k degree of polynomial fit
##' @param iter number of samples to draw from posterior
##' @param lambda set lambda to a specified values
##' @param alpha shape parameter for prior on lambda
##' @param rho rate parameter for prior on lambda
##' @param pb logical indicating weather or not to show a progress bar
##' @aliases btf
##' @author Edward A. Roualdes
##' @export
btf <- function(y='vector', x=NULL, k='int', iter=5e3, alpha=2, rho=0.01, lambda=NULL, pb=TRUE) {

    ## checks
    n <- length(y)
    if (missing(x)) x <- seq_len(n)/n
    nx <- length(x)
    if (nx != n) stop("length of x and y differ.")
    if (any(diff(x) <= 1e-8)) stop("elements of x must be unique")
    D <- genDelta(n, k, x)
    
    ## load module
    loadModule('ind', TRUE)
    individual <- new(ind, y, k, D, alpha, rho) # init new individual

    ## handle lambda
    if (missing(lambda)) {
        f <- function(s) fn(state=s)    # lambda as parameter
    } else {
        if (lambda < 0) {
            print('lambda must be positive; btf will treat lambda as a parameter.')
            f <- function(s) fn(state=s)
        } else {
            f <- function(s) fn(state=s, fixLambda=lambda*lambda) # lambda in model squared
        }
    }

    ## run chains
    chain <- gibbs(iter, f, individual, toCoda, pb=pb)
    attr(chain, 'y') <- y
    attr(chain, 'x') <- x
    class(chain) <- c('btf', class(chain))
    chain
}

fn <- function(state, fixLambda=-1.0) {
    state$upParams(fixLambda)
    c('beta' = state$beta, 's2' = state$s2,
      'l2' = state$l2, 'omega2' = state$o2)
}
