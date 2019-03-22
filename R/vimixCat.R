#' Variational Bayesian inference for unsupervised clustering, mixture of categorical variables
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
#' @export
#'
vimixCatGauss = function(X, K, prior, init = "kmodes", tol = 10e-5,
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
