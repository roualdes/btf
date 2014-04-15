##' Bayesian trend filtering via Eigen
##'
##' @param y response vector
##' @param k degree of polynomial fit
##' @param iters number of samples
##' @author Edward A. Roualdes
##' @export
btf <- function(x=NULL, y='vector' , k='int', iters=5e3, lambda=NULL, rho=0,
                delta=0, l2prior=c('gamma', 'ht'), pb=TRUE) {
    loadModule('ind', TRUE)
    n <- length(y)
    if (missing(x)) {
        x <- seq_len(n)/n
    }
    ## D <- genDelta(n, k, x) # need to re-write this
    D <- genDelta(n, k)
    ## print(head(D[1,]))
    individual <- new(ind, y, k, D, rho, delta) # init new individual
    l2prior <- match.arg(l2prior)               # doesn't yet do anything
    if (missing(lambda)) {
        f <- function(s) fn(state=s)    # lambda as parameter
    } else {
        if (lambda < 0) {
            print('lambda must be positive; fit will treat lambda as a parameter.')
            f <- function(s) fn(state=s)
        } else {
            f <- function(s) fn(state=s, fixLambda=lambda*lambda) # fix lambda
        }
    }
    chain <- gibbs(iters, f, individual, toCoda, pb=pb)
    attr(chain, 'ind') <- individual # output all information possible
    class(chain) <- c('btf', class(chain))
    chain
}

fn <- function(state, fixLambda=-1.0) {
    state$upParams(fixLambda)
    c('beta' = state$beta, 's2' = state$s2,
      'l2' = state$l2, 'omega2' = state$o2)
}
