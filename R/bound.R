
#' Lower bound for multivariate Gaussian
#' @param X NxD matrix data
#' @param model Model parameters
#' @param prior Model parameters
#' @param indep Independent covariates (boolean)
#' @return Value of lower bound
#' @export
boundGauss = function(X, model, prior, indep){

    alpha0 = prior$alpha
    beta0 = prior$beta
    v0 = prior$v
    logW0 = prior$logW
    m0 = prior$m
    M0 = prior$M

    alpha = model$alpha
    beta = model$beta
    v = model$v
    logW = model$logW
    Resp = model$Resp
    m = model$m
    M = model$M

    N = dim(X)[1]
    D = dim(X)[2]
    K = dim(Resp)[2]

    Epz = 0

    Resp[which(Resp<.Machine$double.xmin)] = .Machine$double.xmin # TO DO controllare che questo non alteri i risultati
    Eqz = sum(Resp*log(Resp)) # (10.75) univariate

    logCalpha0 = lgamma(K*alpha0)-K*lgamma(alpha0)
    Eppi = logCalpha0 # (10.73)
    logCalpha = lgamma(sum(alpha))-sum(lgamma(alpha))
    Eqpi = logCalpha # (10.76)

    if(indep){

        Epmu = 0.5*D*K*log(beta0/(2*pi)) - 0.5*sum(beta0/beta) - 0.5*beta0*sum((m-m0)*v/M)# (10.74) 1/2
        Eqmu = 0.5*D*sum(log(beta/(2*pi))) - 0.5*D*K # (10.77) 1/2
        # I have simplified the + 0.5*sum(digamma(v)) + sum(log(1/(0.5*M)))

        vk_wkd = 0; for(k in 1:K) vk_wkd = vk_wkd + sum(v[k]/M[,k])
        EpLambda = (0.5*v0 - 1)*sum(digamma(v)) + sum(log(1/(0.5*M))) -0.5*sum(M0*vk_wkd) # (10.74) 2/2

        EqLambda = 0; for(k in 1:K) EqLambda = EqLambda + sum((0.5*v[k]-1)*sum(digamma(v[k]) + log(1/(2*M[,k])))) -0.5*v[k]# (10.77) 2/2
        # I have simplified the -log(gamma(v0/2))*K*D + 0.5*v0*log(0.5*M0)

        EpX = 0 # (10.71)
        for(n in 1:N){
            for(k in 1:K){
                for(d in 1:D){
                    EpX = EpX + Resp[n,k]*sum(digamma(v[k]) + log(1/(2*M[,k])) - 1/beta[k] - (v[k]*(X[n,d]-m[d,k])^2)/M[d,k] - log(2*pi))
                }
            }
        }
        EpX = 0.5*EpX
    }else{
        Epmu = 0.5*D*K*log(beta0) # (10.74) 1/2
        Eqmu = 0.5*D*sum(log(beta)) # (10.77) 1/2

        logB0 = -0.5*v0*(logW0+D*log(2)); for(d in 1:D) logB0 = logB0 - lgamma(0.5*(v0+1-d))
        EpLambda = K*logB0 # (10.74) 2/2
        logB =  -0.5*v*(logW+D*log(2)); for(k in 1:K) for(d in 1:D) logB[k] = logB[k] - lgamma(0.5*(v[k]+1-d))
        EqLambda = sum(logB) # (10.77) 2/2

        EpX = -0.5*D*N*log(2*pi) # (10.71)
    }

    L = Epz - Eqz + Eppi - Eqpi + Epmu - Eqmu + EpLambda - EqLambda + EpX
    if(!is.finite(L)) stop("Lower bound is not finite")
    L
}

#' Check convergence of variational algorithm
#' @param L Lower bound
#' @param iter Iteration number
#' @param tol Tolerance
#' @param maxiter Maximum number of iterations
#' @param verbose Boolean flag that, if TRUE, prints the outcome of the convergence check
#' @return Boolean flag that indicates if algorithm has converged
#' @export
check_convergence = function(L, iter, tol, maxiter, verbose){
    if (abs(L[iter*2+1]-L[iter*2-1]) < tol*abs(L[iter*2+1])){ # stopping criterion
        if(verbose) message(sprintf("Converged in %d steps.\n", iter))
        conv = TRUE
    }else{
        conv = FALSE
        if(verbose) message(sprintf("Lower bound %f. \n", L[iter*2+1]))
    }
    if (iter == maxiter) warnings(sprintf("Not converged in %d steps.\n", maxiter))
    conv
}
