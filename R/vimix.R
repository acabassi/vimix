#' Variational Bayesian inference for unsupervised clustering
#'
#' @param X NxD data matrix.
#' @param K (Maximum) number of clusters.
#' @param prior Prior parameters (optional).
#' @param init Initialisation method (optional). If it is a vector, it is interpreted as the vector of initial cluster allocations.
#' If it is a string, it is interpreted as the name of the clustering algorithm used for the initialisation (only "kmeans")
#' available at the moment).
#' @param tol Tolerance on lower bound. Default is 10e-20.
#' @param maxiter Maximum number of iterations of the VB algorithm. Default is 2000.
#' @param verbose Boolean flag which, if TRUE, prints the iteration numbers. Default is FALSE.
#' @return A list containing L, the lower bound at each step of the algorithm, label, a vector containing the cluster
#' labels, model, a list containing the trained model structure, and a vector called n_ comp which, if model selection
#' is required, contains the number of mixture components at every step of the VB algorithm.
#' @references Bishop, C.M., 2006. Pattern recognition and machine learning. Springer.
#' @export
#'

vimix = function(X, K, prior, init = NULL, tol = 10e-20,
                 maxiter = 2000, verbose = F){

    if(verbose) message(sprintf("Mixture of multivariate Gaussians \n"))

    if(is.vector(X)){
        N = length(X)
        D = 1
    }else{
        N = dim(X)[1]
        D = dim(X)[2]
    }
    L = rep(-Inf,maxiter*2+1)

    Resp = init(X, K, init, verbose)

    if(missing(prior)){ # set default prior
        if(is.vector(X)){ # if D = 1
            prior = list(alpha = 1, beta = 1, m = mean(X), v = 2, M = 1)
            prior$logW = log(1/prior$M)
        }else{ # if D > 1 and the covariates are not independent
            prior = list(alpha = 1, beta = 1, m = colMeans(X), v = D+1, M = diag(1,D,D))
            prior$logW = log(1/det(prior$M))
        }
    }

    model = prior; model$Resp = Resp # model initialisation
    model = maximizeGauss(X, model, prior)

    for (iter in 1:maxiter){
        if(verbose) message(sprintf("Iteration number %d. ", iter))

        model = expectGauss(X, model) # Expectation step
        L[iter*2] = boundGauss(X, model, prior)/N # Lower bound

        model = maximizeGauss(X, model, prior) # Maximisation step
        L[iter*2+1] = boundGauss(X, model, prior)/N # Lower bound

        if(check_convergence(L, iter, tol, maxiter, verbose)) break # check for convergence
    }
    output = list(L = L[1:iter*2+1], label = apply(model$R, 1, which.max), model=model)

    return(output)
}
