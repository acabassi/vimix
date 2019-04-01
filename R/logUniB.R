#' Compute univariate equivalent of logB
#' @param W Da rivedere
#' @param nu Da rivedere
logUniB <- function(W, nu){
    if(is.na(lgamma(0.5 * nu))) stop("lgamma(0.5*nu) is NA")
    return( 0.5 * nu * log(0.5 / W) - lgamma(0.5 * nu) )
}