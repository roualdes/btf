foobaz <- function(n, k) {
    loadModule('ind', TRUE)
    D <- genDelta(n, k)
    ind <- new(ind, rnorm(10), k, D, 0, 0)
    ## loadModule('individual', TRUE)
    ## ind <- new(individual, rnorm(10), 1, 1, 1) # initialize new individual
    ## Chol
    
    ## eta <- rexp(n-k-1, 0.5)
    ## SigmaInv <- calcSigmaInv(D, eta, n-k-1)
    ## list('CH' = ind$Chol(SigmaInv), 'Sig' = SigmaInv)
    ## Another
    ## ind$Another()
    ## zrnorm
    list('rndGamma' = ind$rndGamma,
         'rInvGauss' = ind$rInvGauss,
         'upParams' = ind$upParams)
}
