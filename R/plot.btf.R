##' plot btf object
##'
##' Plots a btf object with optional credible intervals.
##' @param x btf object
##' @param t domain of function
##' @param burn size of burn-in,
##' @param ptest function specifying how draws from the posterior are summarized
##' @param prob Numeric scalar in (0,1) giving target probability content of the intervals.
##' @param ... extra arguments
##' @author Edward A. Roualdes
##' @aliases plot.btf
##' @export
plot.btf <- function(x, t=NULL, burn = 1e3, ptest=median,
                     prob = 0.95, ...) {
    btf <- x
    y <- attr(btf, 'y')
    n <- length(y)
    if (class(y) == 'ts') {
        t <- as.numeric(time(y))
    } else {
        t <- (1:n)/n
    }

    est <- match.fun(ptest)
    fits <- getPostEst(btf, est=ptest)
    beta <- getPost(btf)
    q <- coda::HPDinterval(beta, prob=prob)
    ry <- range(y)
    rq <- range(q)

    plot(t, y, xlim=range(t), ylim=c(min(ry[1], rq[1]), max(ry[2], rq[2])),
         xlab = 'x', ylab='f(x)', ...)
    lines(t, fits, col="#1AB2FF", lwd=2)
    lines(t, q[,1], col = "#1AB2FF", lty=2, lwd=2)
    lines(t, q[,2], col = "#1AB2FF", lty=2, lwd=2)        
    
}
