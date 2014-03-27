genDelta <- function(n, k) {
    d <- d0 <- bandSparse(n, m=n, c(0,1), diagonals = list(rep(-1, n), rep(1,n-1)))
    if (k > 0) for (i in seq_len(k)) d <- d0 %*% d
    d[seq_len(n-k-1),]
}


matrixClass <- function(M = 'matrix', where = 'string') {
    cl <- class(M)
    print(sprintf("At %s matrix has class %s", where, cl))
}


calcSigmaInv <- function(D, eta, nk) {
    Matrix::t(D) %*% .sparseDiagonal(nk, x=eta) %*% D
}
