genDelta <- function(n, k) {
    d <- d0 <- bandSparse(n, m=n, c(0,1), diagonals = list(rep(-1, n), rep(1,n-1)))
    if (k > 0) for (i in seq_len(k)) d <- d0 %*% d
    d[seq_len(n-k-1),]
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

## colors stolen from
## http://geography.uoregon.edu/datagraphics/color_scales.htm
cols <- c('orange1' = "#FFBF80", 'orange2' = "#FF8000",
          'yellow1' = "#FFFF99", 'yellow2' = "#FFFF33",
          'green1' = "#B2FF8C", 'green2' = "#33FF00",
          'blue1' = "#A6EDFF", 'blue2' = "#1AB2FF",
          'purple1' = "#CCBFFF", 'puple2' = "#664CFF",
          'red1' = "#FF99BF", 'red2' = "#E61A33")


##' get posterior estimates from btf object
##'
##' @param btf btf object
##' @param parameter name of the paramater of interest
##' @param burn number of sample to discard; defaults to 1e3
##' @param est estimate of the object; defaults to median
getPostEst <- function(btf, parameter = c('beta', 's2', 'lambda2', 'omega2'),
                       burn = 1e3, est = median) {
    parameter <- match.arg(parameter)
    est <- match.fun(est)
    idx <- grep(substring(parameter,1,1), varnames(btf))
    M <- dim(btf)[1]
    samples <- as.matrix(btf[(burn+1):M,idx])
    apply(samples,2,est)
}
