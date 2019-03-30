#' Variational Bayesian inference for unsupervised clustering
#'
#' @param X NxD data matrix.
#' @param K (Maximum) number of clusters.
#' @param prior Prior parameters (optional).
#' @param indep Booleand indicator. If TRUE, the features are considered to be independent. Default is FALSE.
#' @param init Initialisation method (optional). If it is a vector, it is interpreted as the vector of initial
#' cluster allocations. If it is a string, it is interpreted as the name of the clustering algorithm used for
#' the initialisation (only "kmeans" and "random") available at the moment).
#' @param select Boolean flag. If TRUE, variable selection is used. Default is FALSE.
#' @param tol Tolerance on lower bound. Default is 10e-20.
#' @param maxiter Maximum number of iterations of the VB algorithm. Default is 2000.
#' @param verbose Boolean flag which, if TRUE, prints the iteration numbers. Default is FALSE.
#' @return A list containing L, the lower bound at each step of the algorithm, label, a vector containing the
#' cluster labels, model, a list containing the trained model structure.
#' @author Alessandra Cabassi \email{alessandra.cabassi@mrc-bsu.cam.ac.uk}
#' @references Bishop, C.M., 2006. Pattern recognition and machine learning. Springer.
#' @examples
#' library(mvtnorm)
#' data <- rbind(rmvnorm(100,c(-3,0)), rmvnorm(100,c(3,0)))
#' output <- vimix(data, 2)
#' @export
#'
vimix = function(X, K, prior, indep = F, init = "kmeans", select = F,
                 tol = 10e-5, maxiter = 2000, verbose = F){

    if(is.vector(X)){
        output = vimixUniGauss(X, K, prior, init, tol, maxiter, verbose)
    }else if(is.list(X) & select){
        stop('Variable selection for categorical data has not been implemented yet.')
    }else if(is.list(X)){
        output = vimixCat(X, K, prior, init, tol, maxiter, verbose)
    }else if(indep & !select){
        output = vimixIndGauss(X, K, prior, init, tol, maxiter, verbose)
    }else if(select){
        output = vimixSelGauss(X, K, prior, init, tol, maxiter, verbose)
    }else{
        output = vimixMulGauss(X, K, prior, init, tol, maxiter, verbose)
    }

    return(output)
}
