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
vimixSelGauss = function(X, K, prior, init = "kmeans", tol = 10e-20,
                         maxiter = 2000, verbose = F){

    if(verbose) message(sprintf("Mixture of univariate Gaussians \n"))

    N = dim(X)[1]
    D = dim(X)[2]

    L = Cl = rep(-Inf, maxiter)

    if(missing(prior)){ # set default prior
            prior = list(alpha = 1/K, beta = 1, m = colMeans(X), v = D + 50, W = rep(1,D), o = 0.5)
    }

    # model initialisation
    Wreshape = as.matrix(prior$W)
    dim(Wreshape) = c(D,1)
    model = list(alpha = rep(prior$alpha,K),
                 beta = matrix(prior$beta, D, K),
                 m =  t(stats::kmeans(X, K, nstart = 25)$centers), # DxK matrix
                 v = matrix(prior$v, D, K),
                 W = Wreshape[, rep(1,K)], # DxK matrix
                 c = rep(.5, D))

    
    lnf = matrix(NA, N, D)
    for(d in 1:D){
        lnf[,d] = log(stats::dnorm(X[,d], mean = prior$m[d], sd = apply(X, 2, stats::sd)[d]))
    }
    model$lnf = lnf

    bound_list = list()
    model_list = list()
    
    for (iter in 2:maxiter){
        if(verbose) message(sprintf("Iteration number %d. ", iter))

        model = expectSelGauss(X, model) # Expectation step
        model = maximizeSelGauss(X, model, prior) # Maximisation step
        model_list[[iter-1]] = model
        bound_list[[iter-1]] = boundSelGauss(X, model, prior) # Lower bound
        L[iter] = bound_list[[iter-1]]$L/N
        Cl[iter] = sum(colSums(model$Resp) > 10e-10*N) # Non-empty clusters

        if(check_convergence(L, iter, tol, maxiter, verbose)) break # check for convergence
    }

    output = list(L = L[1:iter], Cl = Cl[1:iter],
                  label = apply(model$R, 1, which.max), model = model,
                  bound_list = bound_list, model_list = model_list)
}
