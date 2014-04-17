##' plot btf obects
##'
##' @param btf btf object
##' @param x domain of function; defualt is c(0,1)
##' @param burn size of burn-in,
##' @param probs numeric 2-vector of probabilities with values in [0,1]
##' @author Edward A. Roualdes
plot.btf <- function(btf, x = c(0,1), burn = 1e3, probs = c(0.05, 0.95),
                     est = median) {
    y <- attr(btf, 'y')
    n <- length(y)
    x <- seq(x[1], x[2], length.out=n)
    est <- match.fun(est)
    beta <- getPostEst(btf, 'beta', est=est)      # point estimates
    s <- sqrt(getPostEst(btf, 's2', est=est)) # scale estimate
    q <- sapply(beta, function(x, sig) {
        qnorm(probs, x, sig)}, s)
    plot(NULL, xlim = range(x), ylim = range(q))
    plot(y~x, col = 'black', xlab = 'x', ylab='f(x)')
    lines(beta~x, col=cols['blue2'], lwd=2)
    lines(q[1,]~x, col = cols['blue2'], lty=2, lwd=2)
    lines(q[2,]~x, col = cols['blue2'], lty=2, lwd=2)
}
