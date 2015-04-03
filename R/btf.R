##' ##' Bayesian trend filtering via Eigen
##'
##' Fits Bayesian trend filtering hierarchical model to univariate function.
##' Two conditional priors are available: double exponential or generalized double Pareto. 
##'
##' @param y response vector
##' @param x inputs corresponding to y observations
##' @param k degree of polynomial fit
##' @param iter number of samples to draw from posterior
##' @param cond.prior choose the conditional prior on f|sigma
##' @param alpha shape parameter for prior on lambda
##' @param rho rate parameter for prior on lambda
##' @param m sample f every mth iteration, default is m=1 
##' @param D linear transformation of coefficients inside penalty
##' @param debug boolean telling btf to check for NaNs or not
##' @aliases btf
##' @author Edward A. Roualdes
##' @seealso \code{\link[genlasso]{trendfilter}}
##' @references R. J. Tibshirani. Adaptive piecewise polynomial estimation via trend filtering. The Annals of Statistics, 42(1):285-323, 2014.
##' @examples
##' # Cubic trend filtering
##' # from genlasso::trendfilter
##' \dontrun{n <- 100
##' beta0 = numeric(100)
##' beta0[1:40] <- (1:40-20)^3
##' beta0[40:50] <- -60*(40:50-50)^2 + 60*100+20^3
##' beta0[50:70] <- -20*(50:70-50)^2 + 60*100+20^3
##' beta0[70:100] <- -1/6*(70:100-110)^3 + -1/6*40^3 + 6000
##' beta0 <- -beta0
##' beta0 <- (beta0-min(beta0))*10/diff(range(beta0))
##' y <- beta0 + rnorm(n)
##' bfit <- btf(y=y, k=3)
##' plot(bfit, col='grey70')}
##' @export
btf <- function(y='vector', x=NULL, k='int', iter=1e4, cond.prior=c('gdp', 'dexp'), alpha=NULL, rho=NULL, D='Matrix', m=1, debug=FALSE) {

    ## checks
    if (missing(y)) stop('Must provide response vector y.')
    if (!is.numeric(y)) stop('Repsonse vector y must be numeric.')
    if (any(diff(x) == 0)) stop("Elements of x must be unique.")
    if (is.unsorted(x)) stop("X must be in increasing order.")
    n <- length(y)
    if (missing(x)) x <- seq_len(n)/n
    nx <- length(x)
    if (nx != n) stop("Length of x and y differ.")
    if (missing(D) && !missing(k)) {
        if (k < 0 || round(k) != k) stop("Order k must be nonnegative integer.")
        if (k>3) warning(paste("For numerical stability, do not run Bayesian trend filtering with a polynomial order larger than 3."))
        D <- genDelta(n, k, x)
    }else if (missing(k)) {
        stop("Need specify k.")
    }
    nk1 <- nrow(D)
    
    ## which conditional prior?
    cond.prior <- match.arg(cond.prior)
    if ( cond.prior == 'gdp' ) {
        if ( missing(alpha) ) alpha <- 1 else alpha <- alpha  
        if ( missing(rho) ) rho <- 1e-2 else rho <- rho
        
        ## run sampler
        chains <- gdp(iter, y, D, alpha, rho, m, debug)

    } else if ( cond.prior == 'dexp' ) {
        if ( missing(alpha) ) alpha <- 1 else alpha <- alpha
        if ( missing(rho) ) rho <- 1e-4 else rho <- rho

        ## run sampler
        chains <- dexp(iter, y, D, alpha, rho, m, debug)

    } else {
        stop("specified value of cond.prior not understood.")
    }

    ## append some shit for plotting
    names(chains) <- c('f', 's2', 'lambda', 'omega')
    attr(chains, 'y') <- y
    attr(chains, 'x') <- x
    class(chains) <- c('btf', class(chains))
    chains
}
