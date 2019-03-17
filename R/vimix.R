#' Variational Bayesian inference for unsupervised clustering, mixture of categorical variables
#'
#' @param X NxD data matrix.
#' @param K (Maximum) number of clusters.
#' @param prior Prior parameters (optional).
#' @param tol Tolerance on lower bound. Default is 10e-20.
#' @param maxiter Maximum number of iterations of the VB algorithm. Default is 2000.
#' @param verbose Boolean flag which, if TRUE, prints the iteration numbers. Default is FALSE.
#' @return A list containing L, the lower bound at each step of the algorithm, label, a vector containing the
#' cluster labels, model, a list containing the trained model structure.
#' @author Alessandra Cabassi \email{alessandra.cabassi@mrc-bsu.cam.ac.uk}
#' @references Bishop, C.M., 2006. Pattern recognition and machine learning. Springer.
#' @export
#'
vimixCatGauss = function(X, K, prior, init = "kmeans", tol = 10e-5,
                            maxiter = 2000, verbose = F){

    if(verbose) message(sprintf("Mixture of univariate Gaussians \n"))

    if(sum(is.na(X))>0) message("NAs will be treated as additional categories.")

    X <- preprocessCategoricalDataset(X) # convert categories to numbers

    N = dim(X)[1]
    D = dim(X)[2]

    maxNCat <- max(X, na.rm=TRUE) # max number of categories among all features
    nCat <- as.vector(apply(X, 2, max))

    L = Cl = rep(-Inf,maxiter*2)

    if(missing(prior)){ # set default prior
        prior = list(alpha = 1/K)
        prior$eps = matrix(0, D, maxNCat)
        for(d in 1:D){
            prior$eps[d,1:nCat[d]] = 1/nCat[d]
        }
    }

    # model initialisation
    labelInit <- klaR::kmodes(X, modes = K)$cluster
    EPSreshape = prior$eps
    dim(EPSreshape) = c(1,D,maxNCat)
    model = list(alpha = rep(prior$alpha, K),
                 eps = EPSreshape[rep(1,K),,])
    for(i in 1:D){
        for(j in 1:nCat[i]){
            for(k in 1:K){
                model$eps[k,i,j] = prior$eps[i,j] + sum((X[,i]==j)*(labelInit==k))
            }
        }
    }

    for (iter in 2:maxiter){
        if(verbose) message(sprintf("Iteration number %d. ", iter))

        model = expectCatGauss(X, model) # Expectation step
        L[iter*2-1] = boundCatGauss(X, model, prior)/N # Lower bound
        model = maximizeCatGauss(X, model, prior) # Maximisation step
        L[iter*2] = boundCatGauss(X, model, prior)/N # Lower bound
        Cl[iter] = sum(colSums(model$Resp) > 10e-10*N) # Non-empty clusters

        if(check_convergence(L, iter*2, tol, maxiter, verbose)) break # check for convergence
    }

    output = list(L = L[1:(iter*2)], Cl = Cl[1:iter],
                  label = apply(model$R, 1, which.max), model=model)
}

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
vimixIndGauss = function(X, K, prior, init = "kmeans", tol = 10e-20,
                            maxiter = 2000, verbose = F){

    if(verbose) message(sprintf("Mixture of univariate Gaussians \n"))

    N = dim(X)[1]
    D = dim(X)[2]

    L = Cl = rep(-Inf,maxiter)

    if(missing(prior)){ # set default prior
        prior = list(alpha = 1/K, beta = 1, m = mean(X), v = 50 + D, W = rep(100,D))
    }

    # model initialisation
    Wreshape = as.matrix(prior$W)
    dim(Wreshape) = c(D,1)
    model = list(alpha = rep(prior$alpha,K),
                 beta = rep(prior$beta, K),
                 m =  t(stats::kmeans(X, K, nstart = 25)$centers), # DxK matrix
                 v = rep(prior$v, K),
                 W = Wreshape[,rep(1,K)])# DxK matrix

    for (iter in 2:maxiter){
        if(verbose) message(sprintf("Iteration number %d. ", iter))

        model = expectIndGauss(X, model) # Expectation step
        model = maximizeIndGauss(X, model, prior) # Maximisation step
        L[iter] = boundIndGauss(X, model, prior)/N # Lower bound
        Cl[iter] = sum(colSums(model$Resp) > 10e-10*N) # Non-empty clusters

        if(check_convergence(L, iter, tol, maxiter, verbose)) break # check for convergence
    }

    output = list(L = L[1:iter], Cl = Cl[1:iter],
                  label = apply(model$R, 1, which.max), model=model)
}


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

#' Variational Bayesian inference for unsupervised clustering, mixture of multivariate Gaussians
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
#' output <- vimixMulGauss(data, 2)
#' @export
#'
vimixMulGauss = function(X, K, prior, init = "kmeans", tol = 10e-20,
                            maxiter = 2000, verbose = F){

    if(verbose) message(sprintf("Mixture of multivariate Gaussians \n"))

    N = dim(X)[1]
    D = dim(X)[2]

    L = Cl = rep(-Inf,maxiter)

    if(missing(prior)){ # set default prior
        prior = list(alpha = 1/K, beta = 1, m = colMeans(X), v = D+50, W = diag(100,D))
        prior$Winv = solve(prior$W)
    }

    # model initialisation
    Wreshape = prior$W
    dim(Wreshape) = c(dim(Wreshape),1)

    model = list(alpha = rep(prior$alpha,K),
                 beta = rep(prior$beta, K),
                 m =  t(stats::kmeans(X, K, nstart = 25)$centers),
                 v = rep(prior$v, K),
                 W = Wreshape[,,rep(1,K)])

    for (iter in 2:maxiter){
        if(verbose) message(sprintf("Iteration number %d. ", iter))

        model = expectGauss(X, model) # Expectation step
        model = maximizeGauss(X, model, prior) # Maximisation step
        L[iter] = boundGauss(X, model, prior)/N # Lower bound
        Cl[iter] = sum(colSums(model$Resp) > 10e-10*N) # Number of non-empty clusters

        if(check_convergence(L, iter, tol, maxiter, verbose)) break # check for convergence
    }

    output = list(L = L[1:iter], Cl = Cl[1:iter],
                  label = apply(model$R, 1, which.max), model=model)
}

#' Variational Bayesian inference for unsupervised clustering
#'
#' @param X NxD data matrix.
#' @param K (Maximum) number of clusters.
#' @param prior Prior parameters (optional).
#' @param indep Booleand indicator. If TRUE, the features are considered to be independent. Default is FALSE.
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
#' output <- vimix(data, 2)
#' @export
#'
vimix = function(X, K, prior, indep = F, init = "kmeans", tol = 10e-5,
                 maxiter = 2000, verbose = F){

    if(is.vector(X)){
        output = vimixUniGauss(X, K, prior, init, tol, maxiter, verbose)
    }else if(is.list(X)){
        output = vimixCatGauss(X, K, prior, init, tol, maxiter, verbose)
    }else if(indep){
        output = vimixIndGauss(X, K, prior, init, tol, maxiter, verbose)
    }else{
        output = vimixMulGauss(X, K, prior, init, tol, maxiter, verbose)
    }
    return(output)
}
