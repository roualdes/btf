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

