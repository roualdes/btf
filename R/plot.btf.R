##' plot btf obects
##'
##' @author Edward A. Roualdes
plot.btf <- function(x, y, chain = 1, burn = 100, probs = c(0.025, 0.975), ...) {
    betaIdx <- grep('beta', varnames(x))
    n <- length(y)
    t <- seq(0, 1, length.out=n)
    M <- dim(x[,betaIdx])[1]
    samples <- x[(burn+1):M,betaIdx]
    
    plot(1, type='n',
         xlim = c(0,1), ylim = range(y)+c(-1,1),
         xlab = 'x', ylab = 'f(x)')
    points(t, y)
    lines(t, colMeans(samples), col=4)
    
    q <- apply(samples, 2,
               function(x) quantile(x, probs = probs))
    lines(t, q[1,], col='red')
    lines(t, q[2,], col='red')
}
