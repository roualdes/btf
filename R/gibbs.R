##' Gibbs sampler framework
##'
##' implements Gibbs iterations based on updating function fn
##' 
##' @param M number of iterations
##' @param fn function to update parameters in-place
##' @param state member of the population who's state is updated by fn
##' @param toCODA a function to convert state history to CODA object
##' @param pb logical indicating weather or not to show a progress bar
##' @author Edward A. Roualdes
gibbs <- function(M, fn, state, toCODA, pb = TRUE) {
    chain <- vector('list', M)
    if (pb) {
        bar <- txtProgressBar(min=1, max=M, style=3)
        on.exit(close(bar))
    }

    for (m in seq_len(M)) {
        chain[[m]] <- fn(state)
        if (pb) setTxtProgressBar(bar, m)
    }
    toCODA(chain)
}
