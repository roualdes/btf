##' get posterior draws from btf object
##'
##' Return posterior draws from a btf object into a coda object
##' @param btf btf object
##' @param parameter name of the parameter of interest
##' @param burn number of samples to discard
##' @aliases getPost
##' @seealso \code{\link[btf]{getPostEst}}
##' @export
getPost <- function(btf, parameter=c('f', 's2', 'lambda', 'omega'), burn=1e3) {
    parameter <- match.arg(parameter)
    samples <- btf[[parameter]]
    w <- which(samples == 0, arr.ind=TRUE)
    if (length(w)>0) {
        out <- samples[-unique(w[,1]),]
    } else {
        out <- samples
    }
    window(as.mcmc(out), start=burn)
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
getPostEst <- function(btf, parameter=c('f', 's2', 'lambda', 'omega'),
                       burn=1e3, est=NULL) {
    samples <- getPost(btf, parameter=match.arg(parameter), burn=burn)
    if (missing(est)) {
        summary(samples)
    } else {
        apply(as.matrix(samples),2,match.fun(est))
    }
}
