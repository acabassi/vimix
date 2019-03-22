#' Maximization step for univariate Gaussian
#'
#' @param X NxD data matrix
#' @param model Model parameters
#' @param prior Prior parameters
#' @return Updated model parameters
#' @export
#'
maximizeUniGauss = function(X, model, prior){

    alpha0 = prior$alpha
    beta0 = prior$beta
    m0 = prior$m
    v0 = prior$v
    W0 = prior$W
    Resp = model$Resp

    N = length(X)
    D = 1
    K = dim(Resp)[2]

    xbar = m = W = S = rep(0,K)

    Nk = colSums(Resp) + 1e-10 # (10.51)
    alpha = alpha0 + Nk # (10.58)
    beta = beta0 + Nk # (10.60)
    v = v0 + Nk # (10.63)

    for(k in 1:K){
        xbar[k] = sum(Resp[,k]*X)/Nk[k] # (10.52)
        x_cen = X - xbar[k]
        m[k] = (beta0*m0+Nk[k]*xbar[k])/beta[k] # (10.61)
        S[k] = sum((x_cen^2)*Resp[,k])/Nk[k] # (10.53)
        Winv = 1/W0 + Nk[k]*S[k] + ((beta0*Nk[k])/(beta0+Nk[k]))*((xbar[k]-m0)^2) # (10.62)
        W[k] = 1/Winv
    }

    model$alpha = alpha
    model$m = m
    model$W = W
    model$v = v
    model$beta = beta
    model
}
