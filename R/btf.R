##' Bayesian trend filtering via Eigen
##'
##' need better description of this function
##'
##' @param y response vector
##' @param x inputs corresponding to y observations
##' @param k degree of polynomial fit
##' @param iter number of samples to draw from posterior
##' @param cond.prior choose the conditional prior on f|sigma
##' @param alpha shape parameter for prior on lambda
##' @param rho rate parameter for prior on lambda
##' @aliases btf
##' @author Edward A. Roualdes
##' @export
btf <- function(y='vector', x=NULL, k='int', iter=1e4, cond.prior=c('gdP', 'dexp'), alpha=NULL, rho=NULL) {

    ## checks
    n <- length(y)
    if (missing(x)) x <- seq_len(n)/n
    nx <- length(x)
    if (nx != n) stop("length of x and y differ.")
    if (any(diff(x) <= 1e-8)) stop("elements of x must be unique.")
    D <- genDelta(n, k, x)

    ## conditional prior and hyper-parameters
    cond.prior <- match.arg(cond.prior)
    if (tolower(cond.prior) == 'gdp') {
        ## more specific checks
        if (!missing(alpha)) alpha <- alpha else alpha <- 1
        if (!missing(rho)) rho <- rho else rho <- 1
        str <- sprintf('generalized double Pareto conditional prior with alpha = %f and rho = %f',
                       alpha, rho)
        print(str)
        ## run sampler
        stop ('do something with me.')
        
    } else if (tolower(cond.prior) == 'dexp') {

        ## more specific checks
        if (!missing(alpha)) alpha <- alpha else alpha <- 1
        if (!missing(rho)) rho <- rho else rho <- 0
        str <- sprintf('double exponential conditional prior with alpha = %f and rho = %f',
                       alpha, rho)
        print(str)

        ## run sampler
        chain <- dexpBTF(iter, y, k, D, alpha, rho)

        ## tidying
        chain <- as.mcmc(chain)
        varnames(chain) <- c(paste('beta', seq_len(n), sep=''),
                             's2', 'lambda',
                             paste('omega', seq_len(n-k-1), sep=''))

        
    } else {
        stop("specified value of cond.prior not understood.")
    }

    ## append some shit for printing
    attr(chain, 'y') <- y
    attr(chain, 'x') <- x
    class(chain) <- c('btf', class(chain))
    chain
}
