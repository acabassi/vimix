#' Expectation step for mixture of multivariate Gaussians
#'
#' @param X NxD data matrix.
#' @param model Model parameters.
#' @return Updated model parameters.
#' @export
#'
expectGauss = function(X, model){

    alpha = model$alpha

    v = model$v
    beta = model$beta
    m = model$m
    M = model$M

    if(is.vector(X)){
        N = length(X)
        D = 1
    }else{
        N = dim(X)[1]
        D = dim(X)[2]
    }

    K = length(v)

    logRho = matrix(NA, N, K)

    logW = rep(NA,K)
    for (k in 1:K){
        ElnPi = expeCat(alpha, k)
        logW[k] = update_logW(M[,,k])
        logRho[,k] = expeCat(alpha, k) + expeGauss(X, D, N, v[k], m[,k], M[,,k], beta[k], logW[k])
    }
    model$logW = logW

    Resp = logRho_to_Resp(logRho)

    model$Resp = Resp

    model
}

#' Update log(W)
#' @param M = 1/W
#' @export
#'
update_logW = function(M) {
    if(is.matrix(M)){
        logW = log(1/det(M))
    }else{
        logW = log(1/M)
    }
    return(logW)
}

#' Update log(Rho) for independent Gaussians
#' @param X NxD data matrix.
#' @param D Number of covariates.
#' @param N Number of observations.
#' @param v_k Precision hyperparameter.
#' @param m_k Mean hyperparameter.
#' @param M_k Precision hyperparameter.
#' @param beta_k Mean hyperparameter.
#' @param logW_k log(1/det(M)), not needed if indep is TRUE).
#' @export
#'
expeGauss = function(X, D, N, v_k, m_k, M_k, beta_k, logW_k){

    if(D==1){
        # v = model$v_y[k]; m = model$m_y[k]; M = model$M_y[k]; beta = model$beta_y[k]; ms = FALSE
        ElnLa = digamma(0.5*v_k) + log(0.5*M_k) # (10.65) univariate case
        diff = X - m_k
        ExmuLaxmu = 1/beta_k + v_k*(diff^2)/M_k
        logRho = 0.5*ElnLa - 0.5*ExmuLaxmu
    }else{
        # logW_k = logW[k]; v_k = v[k]; m_k = m[,k]; M_k = M[,,k]; beta_k = beta[k]
        ElnLa = D*log(2) + logW_k + digamma(0.5*(v_k+1-c(1:D))) # (10.65)
        diff = sweep(X, 2, m_k, FUN="-")
        xmuLaxmu = rep(NA, N); for(n in 1:N) xmuLaxmu[n] = t(diff[n,])%*%solve(M_k,diff[n,]) # TODO sistemare
        ExmuLaxmu = D/beta_k + v_k*xmuLaxmu # (10.64)
        logRho = 0.5*ElnLa - 0.5*ExmuLaxmu # (10.46)
    }
    return(logRho)
}

#' Expectation of log(pi)
#' @param alpha Dirichlet hyperparameter.
#' @param k Mixture component of interest.
#' @export
expeCat = function(alpha, k){
    ElnPi = digamma(alpha[k]) - digamma(sum(alpha)) # (10.66)
    ElnPi
}


#' Transform log(Rho) matrix to responsibilities matrix
#' @param logRho log(Rho).
#' @return Matrix of responsibilities.
#' @export
logRho_to_Resp = function(logRho){

    maxLogRho = apply(logRho,1,max)
    logRho = logRho - maxLogRho
    logSumExpLogRho =  log(apply(exp(logRho),1,sum))
    logResp = logRho - logSumExpLogRho
    Resp = exp(logResp)
}
