##' Bayesian trend filtering
##'
##' @param y response vector
##' @param k degree of polynomial fit
##' @param iters number of samples
##' @author Edward A. Roualdes
##' @export
tf <- function(y='vector' , k='int', iters=1e4, lambda = NULL,
               l2prior = c('gamma', 'ht')) {
    loadModule('individual', TRUE)
    ind <- new(individual, y, k, 1, 1) # initialize new individual
    l2prior <- match.arg(l2prior)
    upParams <- function(indiv) {
        ## update in place
        indiv$upBeta(); indiv$upOmega();
        if (is.null(lambda)) {
            if (l2prior == 'ht') {
                xstar <- prop()
                a <- 
                alpha <- min(1, a)
            } else {
                indiv$upLambda();                
            }
        } else {
            ind$l2 <- lambda*lambda
        }
        indiv$upSig()
        ## only return what we want a history of
        c('beta' = indiv$beta, 'l2' = indiv$l2,
          'omega2' = indiv$o2, 's2' = indiv$s2)
    }
    chain <- gibbs(iters, upParams, ind, toCoda)
    attr(chain, 'ind') <- individual # output all information possible
    class(chain) <- c('btf', class(chain))
    chain
}

##' update state
##'
##' @param indiv an individual (c++ object) for Bayesian trend filtering
##' @author Edward A. Roualdes
upParams <- function(indiv) {
    ## update in place
    indiv$upBeta()
    indiv$upOmega()
    indiv$upLambda()
    indiv$upSig()
    ## only return what we want a history of
    c('beta' = indiv$beta, 'l2' = indiv$l2,
      'omega2' = indiv$o2, 's2' = indiv$s2)
}

##' transform state's history into CODA object
##'
##' @param chain list of states through history
##' @author Edward A. Roualdes
toCoda <- function(chain) {
    names <- names(chain[[1]])
    hist <- mcmc(t(simplify2array(chain)))
    varnames(hist) <- names
    hist
}
