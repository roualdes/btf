##' bayesian trend filtering
##'
##' @author Edward A. Roualdes

btfR <- function(y, k, iters=1e4, bar=TRUE, rho=NULL, delta=NULL) {

    ## initialize
    n <- length(y)
    nk <- n-k-1                         
    I <- .sparseDiagonal(n)             # identity
    D <- genDelta(n, k)                 # discrete derivative matrix
    l <- n+nk+2                         # length of parameters
    chain <- matrix(0, iters+1, l)      # chain history

    if (missing(rho) ||  missing(delta)) {
        l2 <- rgamma(1, nk, 1)
    } else {
        l2 <- rgamma(1, nk+rho, 1+delta)
    }
    o2 <- rexp(nk, 0.5); eta <- 1/o2    # we model eta = 1/o2
    SigmaInv <- calcSigmaInv(D, eta, nk)
    s2 <- rinvGamma(1, 0.1, 0.1)
    beta <- y

    ## iterate
    iters <- iters+1L
    chain[1,] <- c(beta, s2, l2, o2)
    i <- 1L
    if (bar) {
        pb <- txtProgressBar(i, iters+1, style=3)
        on.exit(close(pb))
    } 

    while ( (i <- i+1) <= iters ) {
        ## update parameters
        beta <- upBeta(y, I, n, SigmaInv, s2)               # beta
        l2 <- rgamma(1, nk, rate=sum(o2)/2)                 # lambda^2
        eta <- rinvGauss(nk, sqrt(l2*s2)/abs(D%*%beta), l2) # eta => 1/o2
        tmp <- as.numeric(crossprod(y-beta) + (t(beta)%*%SigmaInv%*%beta))
        s2 <- rinvGamma(1, n, 2/tmp) # sigma^2

        ## calculate some things before we store chains
        o2 <- 1/eta                          # omega^2
        SigmaInv <- calcSigmaInv(D, eta, nk) # calculate Sigma_f^{-1}

        ## store chains
        chain[i,] <- c(beta, s2, l2, o2)

        if (bar) setTxtProgressBar(pb,i)
    }
    class(chain) <- c('btf', class(chain))
    attr(chain, 'SigmaInv') <- SigmaInv
    attr(chain, 'Delta') <- D
    chain
}
