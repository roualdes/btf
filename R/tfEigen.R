##' Bayesian trend filtering via Eigen
##'
##' @param y response vector
##' @param k degree of polynomial fit
##' @param iters number of samples
##' @author Edward A. Roualdes
##' @export
btfE<- function(y='vector' , k='int', iters=1e4, lambda = NULL, rho = 0,
                delta = 0, l2prior = c('gamma', 'ht')) {
    loadModule('ind', TRUE)
    D <- genDelta(length(y), k)
    individual <- new(ind, y, k, D, rho, delta) # init new individual
    l2prior <- match.arg(l2prior)
    chain <- gibbs(iters, fn, individual, toCoda)
    attr(chain, 'ind') <- individual # output all information possible
    class(chain) <- c('btf', class(chain))
    chain
}

fn <- function(state) {
    state$upParams()
    c('beta' = state$beta, 's2' = state$s2,
      'l2' = state$l2, 'omega2' = state$o2)
}
