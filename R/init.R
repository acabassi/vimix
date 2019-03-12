#' Initialize responisibility matrix
#'
#' This function allows you to do variational Bayesian inference for Gaussian mixture.
#' @param X NxD data matrix
#' @param K Number of mixture components
#' @param init Initialisation method (optional). If it is a vector, it is interpreted as the vector of initial cluster allocations.
#' If it is a string, it is interpreted as the name of the clustering algorithm used for the initialisation (only "kmeans")
#' available at the moment).
#' @param verbose Boolean flag which, if TRUE, prints the type of initialization
#' @return Resp NxK matrix of responsibilities
#' @export
#'
init = function(X, K, init = F, verbose = T){

    if(is.vector(X)){
        N = length(X)
    }else{
        N = dim(X)[1]
    }

    if(!is.vector(init)){
        if(verbose) message(sprintf("Initialization of responsibilities: user-defined labels. \n"))
        label = init # assign some user-defined labels
    }else if(init=="kmeans"){
        if(verbose) message(sprintf("Initialization of responsibilities: k-means. \n"))
        km = stats::kmeans(X, K, iter.max = 100)
        label = km$cluster # initial labels correspond to the output of the k-means algorithm
    }else{
        if(verbose) message(sprintf("Initialization of responsibilities: random assignment to mixture components.\n"))
        label = sample(1:K, N, replace=T) # assign random labels
    }
    Resp = as.matrix(Matrix::sparseMatrix(i = c(1:N), j = label, x = 1)) # create responsibility matrix
    Resp
}
