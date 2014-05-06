##' Bayesian trend filtering via Eigen
##'
##' need better description of this function
##'
##' @param y response vector
##' @param x inputs corresponding to y observations
##' @param k degree of polynomial fit
##' @param iter number of samples to draw from posterior
##' @param cond.prior choose the conditional prior on f|sigma
##' @param lambda set lambda to a specified values
##' @param alpha shape parameter for prior on lambda
##' @param rho rate parameter for prior on lambda
##' @param pb logical indicating weather or not to show a progress bar
##' @aliases btf
##' @author Edward A. Roualdes
##' @export
btf <- function(y='vector', x=NULL, k='int', iter=5e3, cond.prior=c('gdp', 'dexp'), alpha=NULL, rho=NULL, lambda=NULL, pb=FALSE) {

    ## checks
    n <- length(y)
    if (missing(x)) x <- seq_len(n)/n
    nx <- length(x)
    if (nx != n) stop("length of x and y differ.")
    if (any(diff(x) <= 1e-8)) stop("elements of x must be unique.")
    D <- genDelta(n, k, x)

    ## conditional prior and hyper-parameters
    cond.prior <- match.arg(cond.prior)
    if (cond.prior == 'gdp') {
        if (!missing(alpha)) alpha <- alpha else alpha <- 1
        if (!missing(rho)) rho <- rho else rho <- 1
        cprior <- 1
    } else if (cond.prior == 'dexp') {
        if (!missing(alpha)) alpha <- alpha else alpha <- 1
        if (!missing(rho)) rho <- rho else rho <- 0
        cprior <- 2
    } else stop("specified value of cond.prior not understood.")
    print(paste('alpha = ', alpha))
    print(paste('rho = ', rho))
    print(cprior)
    
    ## load module
    loadModule('ind', TRUE)
    individual <- new(ind, y, k, D, cprior, alpha, rho) # init new individual

    ## handle lambda
    if (missing(lambda)) {
        f <- function(s) upParams(state=s)    # lambda as parameter
    } else {
        if (lambda < 0) {
            print('lambda must be positive; btf will treat lambda as a parameter.')
            f <- function(s) upParams(state=s)
        } else {
            f <- function(s) upParams(state=s, fixLambda=lambda*lambda) # lambda in model squared
        }
    }

    ## run chains
    chain <- gibbs(iter, f, individual, toCoda, pb=pb)
    attr(chain, 'y') <- y
    attr(chain, 'x') <- x
    class(chain) <- c('btf', class(chain))
    chain
}

upParams <- function(state, fixLambda=-1.0) {
    state$upParams(fixLambda)
    c('beta' = state$beta, 's2' = state$s2,
      'lambda' = state$l, 'lambda2' = state$l2, 'omega' = state$o)
}
