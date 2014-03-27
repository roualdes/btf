##' plot btf obects
##'
##' @param x coda object that contains a btf history
##' @param y response vector for btf model
##' @param chain number of chains included in x
##' @param burn size of burn-in,
##' @param probs numeric 2-vector of probabilities with values in [0,1]
##' @author Edward A. Roualdes
plot.btf <- function(x, y, chain = 1, burn = 1e3, probs = c(0.025, 0.975),
                     est = c('median', 'mean')) {
    betaIdx <- grep('beta', varnames(x))
    n <- length(y)
    t <- seq(0, 1, length.out=n)
    M <- dim(x[,betaIdx])[1]
    samples <- x[(burn+1):M,betaIdx]

    ## colors stolen from
    ## http://geography.uoregon.edu/datagraphics/color_scales.htm
    cols <- c('orange1' = "#FFBF80", 'orange2' = "#FF8000",
              'yellow1' = "#FFFF99", 'yellow2' = "#FFFF33",
              'green1' = "#B2FF8C", 'green2' = "#33FF00",
              'blue1' = "#A6EDFF", 'blue2' = "#1AB2FF",
              'purple1' = "#CCBFFF", 'puple2' = "#664CFF",
              'red1' = "#FF99BF", 'red2' = "#E61A33")

    q <- apply(samples, 2,
               function(x) quantile(x, probs = probs))
    e <- apply(samples, 2, est)

    plot(y~t, col = 'black', xlab = 'x', ylab='f(x)')
    lines(e~t, col=cols['blue2'], lwd=2)
    lines(q[1,]~t, col = cols['orange2'], lty=2, lwd=2)
    lines(q[2,]~t, col = cols['orange2'], lty=2, lwd=2)
}
