#' Maximization step for multivariate Gaussian
#'
#' @param X NxD data matrix
#' @param model Model parameters
#' @param prior Prior parameters
#' @return Updated model parameters
#' @export
#'

maximizeGauss = function(X, model, prior){

    alpha0 = prior$alpha
    beta0 = prior$beta
    m0 = prior$m
    v0 = prior$v
    M0 = prior$M
    logW0 = prior$logW
    M0 = prior$M
    Resp = model$Resp

    if(is.vector(X)){
        N = length(X)
        D = 1
    }else{
        N = dim(X)[1]
        D = dim(X)[2]
    }
    K = dim(Resp)[2]

    Nk = colSums(Resp) # (10.51)
    beta = beta0 + Nk # (10.60)
    v = v0 + Nk  # (10.63)

    xbar = matrix(0,D,K)
    m = matrix(NA,D,K)

    M = array(NA, c(D,D,K))
    for(k in 1:K){
        xbar[,k] = update_xbar(X, Resp[,k]) # (10.52)
        M[,,k] = update_invW(X, Resp[,k], M0, beta0, m0, xbar[,k]) # (10.62)
        m[,k] = update_m(beta0, m0, Nk[k], xbar[,k], beta[k]) # (10.61)
    }

    alpha0 = prior$alpha
    alpha = alpha0 + Nk # (10.58)
    model$alpha = alpha

    model$m = m
    model$M = M
    model$v = v
    model$beta = beta
    model
}

#' Update xbar
#'
#' @param X NxD data matrix.
#' @param Resp Vector of responsibilities for the mixture component of interest.
#' @export
#'
update_xbar = function (X, Resp){

    Nk = sum(Resp)
    if(is.vector(X)){
        xbarNk = sum(Resp*X)
        xbar = xbarNk/Nk
    }else{
        xbarNk = colSums(sweep(X, MARGIN=1, Resp, `*`))
        xbar = xbarNk/Nk
    }
    if(is.nan(sum(xbar))) stop()
    xbar
}

#' Update M = inv(W)
#'
#' @param X NxD data matrix.
#' @param Resp_k Vector of responsibilities for mixture component k.
#' @param M0_k Prior precision for mixture component k.
#' @param beta0_k Prior beta for mixture component k.
#' @param m0_k Prior m for mixture component k.
#' @param xbar_k Weighted mean for mixture component k.
#' @export
#'
update_invW = function(X, Resp_k, M0_k, beta0_k, m0_k, xbar_k){

    # Resp_k = Resp[,k]; M0_k = M0; beta0_k = beta0; m0_k = m0; xbar_k = xbar[,k]

    Nk = sum(Resp_k)

    if(is.vector(X)){ # Univariate Gaussian
        N = length(X)
        SkNk = 0
        for(n in 1:N) SkNk = SkNk + Resp_k[n]*(X[n]-xbar_k)^2 # (10.53) univariate
        M_k = M0_k + SkNk + (beta0_k*Nk)/(beta0_k+Nk)*(xbar_k-m0_k)^2
    }else{ # Multivariate Gaussian
        N = dim(X)[1]
        D = dim(X)[2]
        SkNk = matrix(0,D,D)
        for(n in 1:N) SkNk = SkNk + Resp_k[n]*tcrossprod(X[n,]-xbar_k) # (10.53)
        M_k = M0_k + SkNk + (beta0_k*Nk)/(beta0_k+Nk)*tcrossprod(xbar_k-m0_k)

    }
    M_k
}

#' Update mean vector
#'
#' @param beta0 Prior beta.
#' @param m0 Prior m.
#' @param Nk Number of observations in component k.
#' @param xbar Weighted mean of component k.
#' @param beta Parameter for precision of mu.
#' @export
#'
update_m = function(beta0, m0, Nk, xbar, beta){
    m = (beta0*m0 + Nk*xbar)/beta
    if(is.nan(sum(m))) stop()
    m
}
