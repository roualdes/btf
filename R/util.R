gendelta <- function(n, k) {
    d <- d0 <- bandSparse(n, m=n, c(0,1), diagonals = list(rep(-1, n), rep(1,n-1)))
    if (k > 0) for (i in seq_len(k)) d <- d0 %*% d
    d[seq_len(n-k-1),]
}


##' generate Matrix Delta^{k+1}
##'
##' @param n sample size
##' @param k order of fit
##' @param x vector of inputs on domain
##' @author Edward A. Roualdes
genDelta <- function(n, k, x) {
    nk <- n-k-1
    d <- Matrix(0, nk, n)
    for (i in seq_len(nk)) {
        ik <- i+k+1
        idx <- i:ik
        tmp <- x[ik] - x[i]
        z <- x[idx]
        d[i,idx] <- tmp*sapply(seq_along(z),function(j) 1/prod(z[j]-z[-j]))
    }
    gamma(k+1)*d/n^k
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

##' extract posterior samples from btf object
##'
##' @param btf object
##' @param parameter name of the parameter of interest
##' @param burn number of samples to discard
##' @export
getPost <- function(btf, parameter = c('beta', 's2', 'lambda', 'omega'),
                    burn = 1e3) {
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
getPostEst <- function(btf, parameter = c('beta', 's2', 'lambda', 'omega'),
                       burn = 1e3, est = median) {
    samples <- getPost(btf, parameter=match.arg(parameter), burn=burn)
    apply(as.matrix(samples),2,match.fun(est))
}
