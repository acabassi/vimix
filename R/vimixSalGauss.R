#' Variational Bayesian inference for unsupervised clustering, mixture of independent Gaussians
#'
#' @param X NxD data matrix.
#' @param K (Maximum) number of clusters.
#' @param prior Prior parameters (optional).
#' @param init Initialisation method (optional). If it is a vector, it is interpreted as the vector of initial
#' cluster allocations. If it is a string, it is interpreted as the name of the clustering algorithm used for
#' the initialisation (only "kmeans" and "random") available at the moment).
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
#' output <- vimixIndGauss(data, 2)
#' @export
#'
vimixSalGauss = function(X, K, prior, init = "kmeans", tol = 10e-20,
                         maxiter = 2000, verbose = F){

    if(verbose) message(sprintf("Mixture of univariate Gaussians \n"))

    N = dim(X)[1]
    D = dim(X)[2]

    L = Cl = rep(-Inf, maxiter)

    if(missing(prior)){ # set default prior
            prior = list(alpha = 1e-16, beta = 1e-16, m = colMeans(X), c = 1e-16)
    }
    
    model = list(alpha = matrix(prior$alpha, D,K),
                 beta = matrix(prior$beta, D, K),
                 m =  t(stats::kmeans(X, K, nstart = 25)$centers), # DxK matrix
                 c = matrix(prior$c, D, K),
                 pi = rep((1/K), K),
                 rho = matrix(.5, N, D), 
                 Resp = matrix(NA, N, K))
    
    for (iter in 2:maxiter){
        if(verbose) message(sprintf("Iteration number %d. ", iter))

        model = updateSalGauss(X, model, prior)
        Cl[iter] = sum(colSums(model$Resp) > 10e-10*N) # Non-empty clusters

    }

    output = list(Cl = Cl[1:iter], label = apply(model$Resp, 1, which.max), model = model)
}
