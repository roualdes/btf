##' extract posterior samples from btf object
##'
##' @param btf object
##' @param parameter name of the parameter of interest
##' @param burn number of samples to discard
##' @export
getPost <- function(btf, parameter=c('beta', 's2', 'lambda', 'omega'),
                    burn=1e3) {
    parameter <- match.arg(parameter)
    idx <- grep(parameter, varnames(btf))
    window(btf[,idx], start=burn)
}

##' get posterior estimates from btf object
##'
##' @param btf btf object
##' @param parameter name of the paramater of interest
##' @param burn number of samples to discard
##' @param est estimate of the object
##' @export
getPostEst <- function(btf, parameter=c('beta', 's2', 'lambda', 'omega'),
                       burn=1e3, est=mean) {
    samples <- getPost(btf, parameter=match.arg(parameter), burn=burn)
    apply(as.matrix(samples),2,match.fun(est))
}
