##' get posterior draws from btf object
##'
##' Return posterior draws from a btf object into a coda object
##' @param btf btf object
##' @param parameter name of the parameter of interest
##' @param burn number of samples to discard
##' @aliases getPost
##' @seealso \code{\link[btf]{getPostEst}}
##' @export
getPost <- function(btf, parameter=c('beta', 's2', 'lambda', 'omega', 'alpha', 'rho'),
                    burn=1e3) {
    parameter <- match.arg(parameter)
    idx <- grep(parameter, varnames(btf))
    window(btf[,idx], start=burn)
}

##' get posterior estimates from btf object
##'
##' Return posterior estimates, computed from btf posterior draws.
##' User can supply their own estimating function, \code{est}, applied
##' to the columns of \code{\link[btf]{getPost}}.
##' @param btf btf object
##' @param parameter name of the paramater of interest
##' @param burn number of samples to discard
##' @param est function to estimate from the posterior samples of \code{parameter}
##' @aliases getPostEst
##' @seealso \code{\link[btf]{getPost}}
##' @export
getPostEst <- function(btf, parameter=c('beta', 's2', 'lambda', 'omega', 'alpha', 'rho'),
                       burn=1e3, est=NULL) {
    samples <- getPost(btf, parameter=match.arg(parameter), burn=burn)
    if (missing(est)) {
        if (is.matrix(samples)) {
            samps <- as.mcmc.list(samples)
        } else {
            samps <- as.mcmc(samples)
        }
        summary(samps)
    } else {
        apply(as.matrix(samples),2,match.fun(est))
    }
}
