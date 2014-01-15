##' Gibbs sampler framework
##'
##' @param M number of iterations
##' @param fn function to update parameters in-place
##' @param individual member of the population who's state is updated by fn
##' @param toCODA a function to convert state history to CODA object
##' @author Edward A. Roualdes
gibbs <- function(M, fn, individual, toCODA) {
    chain <- vector('list', M)
    pb <- txtProgressBar(min=1, max=M, style=3); on.exit(close(pb))
    for (m in seq_len(M)) {
        chain[[m]] <- fn(individual)
        setTxtProgressBar(pb, m)
    }
    toCODA(chain)
}
