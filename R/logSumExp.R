#' Use the log sum exp trick to improve numerical stability
#' @param x Vector of values
log_sum_exp <- function(x) {
    # Computes log(sum(exp(x))
    offset <- max(x)
    s <- log(sum(exp(x - offset))) + offset
    i <- which(!is.finite(s))
    if (length(i) > 0) { s[i] <- offset }
    return(s)
}
