upBeta <- function(obs, eye, n, siginv, sig2) {
    CH <- Cholesky(as(eye+siginv, 'symmetricMatrix'))
    as.vector(zrnormChol(n, as.vector(solve(CH, obs)), sig2*solve(CH)))
}
