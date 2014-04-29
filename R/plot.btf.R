##' plot btf obects
##'
##' @param x btf object
##' @param t domain of function
##' @param burn size of burn-in,
##' @param probs numeric 2-vector of probabilities with values in [0,1]
##' @param est function specifying how draws from the posterior are summarized
##' @param ... extra arguments to be passed as methods
##' @author Edward A. Roualdes
##' @aliases plot.btf
##' @export
plot.btf <- function(x, t=NULL, burn = 1e3, probs = c(0.05, 0.95),
                     est = median, ...) {
    btf <- x
    y <- attr(btf, 'y')
    n <- length(y)
        if (missing(t)) {
        t <- (1:n)/n
    }

    est <- match.fun(est)
    beta <- getPostEst(btf, 'beta', est=est)      # point estimates
    q <- getPostEst(btf, 'beta', est=function(z) quantile(z, c(0.025, 0.975))) # ci
    
    plot(y~t, xlim = range(t), ylim = range(q),
         col = 'black', xlab = 'x', ylab='f(x)', ...)
    lines(beta~t, col=cols['blue2'], lwd=2)
    lines(q[1,]~t, col = cols['blue2'], lty=2, lwd=2)
    lines(q[2,]~t, col = cols['blue2'], lty=2, lwd=2)
}
