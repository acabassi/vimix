#' Use the log sum exp trick to improve numerical stability
#' @param x Vector of values
log_sum_exp <- function(x) {
    # Computes log(sum(exp(x))
    offset <- max(x)
    ss <- log(sum(exp(x - offset))) + offset
    ii <- which(!is.finite(ss))
    if (length(ii) > 0) { ss[ii] <- offset }
    return(ss)
}
