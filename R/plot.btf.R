##' plot btf obects
##'
##' @param x coda object that contains a btf history
##' @param y response vector for btf model
##' @param chain number of chains included in x
##' @param burn size of burn-in,
##' @param probs numeric 2-vector of probabilities with values in [0,1]
##' @author Edward A. Roualdes
plot.btf <- function(x, chain = 1, burn = 1e3, probs = c(0.025, 0.975),
                     est = median) {
    betaIdx <- grep('beta', varnames(x))
    y <- attr(x, 'ind')$y
    n <- length(y)
    t <- seq(0, 1, length.out=n)
    M <- dim(x[,betaIdx])[1]
    est <- match.fun(est)
    samples <- x[(burn+1):M,betaIdx]

    ## replace with getPostEst fn; see TODO:util.R
    q <- apply(samples, 2,
               function(x) quantile(x, probs = probs))
    e <- apply(samples, 2, est)

    plot(y~t, col = 'black', xlab = 'x', ylab='f(x)')
    lines(e~t, col=cols['blue2'], lwd=2)
    lines(q[1,]~t, col = cols['orange2'], lty=2, lwd=2)
    lines(q[2,]~t, col = cols['orange2'], lty=2, lwd=2)
}
