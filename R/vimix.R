#' Variational Bayesian inference for unsupervised clustering
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
#' cluster labels, model, a list containing the trained model structure, and a vector called n_ comp which, if
#' model selection is required, contains the number of mixture components at every step of the VB algorithm.
#' @author Alessandra Cabassi \email{alessandra.cabassi@mrc-bsu.cam.ac.uk}
#' @references Bishop, C.M., 2006. Pattern recognition and machine learning. Springer.
#' @examples
#' ## Load a dataset containing 200 2-dimensional data points
#' data <- as.matrix(read.csv(system.file("extdata", "example1-data.csv", package = "vimix"),
#' row.names = 1))
#'
#' ## Use variational inference for mixture of Gaussians to find clusters
#' output <- vimix(data, 2)
#' @export
#'

vimix = function(X, K, prior, init = "kmeans", tol = 10e-20,
                 maxiter = 2000, verbose = F){

    if(verbose) message(sprintf("Mixture of multivariate Gaussians \n"))

    if(is.vector(X)){
        N = length(X)
        D = 1
    }else{
        N = dim(X)[1]
        D = dim(X)[2]
    }
    L = rep(-Inf,maxiter)

    # Resp = init(X, K, init, verbose)

    if(missing(prior)){ # set default prior
            prior = list(alpha = 1/K, beta = 1, m = colMeans(X), v = D+50, W = diag(100,D))
            prior$Winv = solve(prior$W)
    }

    # model initialisation

    Wreshape = prior$W
    dim(Wreshape) = c(dim(Wreshape),1)

    model = list(alpha = rep(prior$alpha,K),
                 beta = rep(prior$beta, K),
                 m =  t(kmeans(X, K, nstart = 25)$centers),
                 v = rep(prior$v, K),
                 W = Wreshape[,,rep(1,K)])
    # model = maximizeGauss(X, model, prior)

    # cat('model$alpha', model$alpha, '\n')
    # cat('model$beta', model$beta, '\n')
    # cat('model$m', model$m, '\n')
    # cat('model$v', model$v, '\n')
    # cat('model$W', model$W[,1,1], '\n')
    # cat('model$Resp[1,', (model$Resp), '\n')

    for (iter in 2:maxiter){
        if(verbose) message(sprintf("Iteration number %d. ", iter))

        model = expectGauss(X, model) # Expectation step
        # L[iter*2] = boundGauss(X, model, prior)/N # Lower bound

        # cat('model$Resp', head(model$Resp), '\n')

        model = maximizeGauss(X, model, prior) # Maximisation step
        L[iter] = boundGauss(X, model, prior)/N # Lower bound

        # cat('model$alpha', model$alpha, '\n')
        # cat('model$beta', model$beta, '\n')
        # cat('model$m', model$m[,1], '\n')
        # cat('model$v', model$v[1], '\n')
        # cat('model$M', model$W[,1,1], '\n')

        if(check_convergence(L, iter, tol, maxiter, verbose)) break # check for convergence
    }
    output = list(L = L[1:iter], label = apply(model$R, 1, which.max), model=model)

    return(output)
}
