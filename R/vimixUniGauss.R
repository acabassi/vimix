#' Variational Bayesian inference for unsupervised clustering, mixture of univariate Gaussians
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
#' data <- c(rnorm(100, -4), rnorm(100, 4))
#' output <- vimixUniGauss(data, 3)
#' @export
#'
vimixUniGauss = function(X, K, prior, init = "kmeans", tol = 10e-20,
                         maxiter = 2000, verbose = F){

    if(verbose) message(sprintf("Mixture of univariate Gaussians \n"))

    N = length(X)
    D = 1

    L = Cl = rep(-Inf,maxiter)

    if(missing(prior)){ # set default prior
        prior = list(alpha = 1/K, beta = 1, m = mean(X), v = D + 50, W = 100)
    }

    # model initialisation
    model = list(alpha = rep(prior$alpha,K),
                 beta = rep(prior$beta, K),
                 m =  as.vector(stats::kmeans(X, K, nstart = 25)$centers),
                 v = rep(prior$v, K),
                 W = rep(prior$W, K))

    for (iter in 2:maxiter){
        if(verbose) message(sprintf("Iteration number %d. ", iter))

        model = expectUniGauss(X, model) # Expectation step
        model = maximizeUniGauss(X, model, prior) # Maximisation step
        L[iter] = boundUniGauss(X, model, prior)/N # Lower bound
        Cl[iter] = sum(colSums(model$Resp) > 10e-10*N)

        if(check_convergence(L, iter, tol, maxiter, verbose)) break # check for convergence
    }

    output = list(L = L[1:iter], Cl = Cl[1:iter],
                  label = apply(model$R, 1, which.max), model=model)
}
