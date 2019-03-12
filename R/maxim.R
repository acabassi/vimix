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
    W0 = prior$W
    W0inv = prior$Winv
    Resp = model$Resp

    N = dim(X)[1]
    D = dim(X)[2]
    K = dim(Resp)[2]

    xbar = m = matrix(0,D,K)
    W = S = array(NA, c(D,D,K))

    Nk = colSums(Resp) + 1e-10 # (10.51)
    alpha = alpha0 + Nk # (10.58)
    beta = beta0 + Nk # (10.60)
    v = v0 + Nk # (10.63)

    for(k in 1:K){
        xbar[,k] = (Resp[,k]%*%X)/Nk[k] # (10.52)
        x_cen = sweep(X, MARGIN = 2, STATS = xbar[,k], FUN = "-")
        S[,,k] = t(x_cen)%*%(x_cen*Resp[,k])/Nk[k] # (10.53)
        m[,k] = (beta0*m0+Nk[k]*xbar[,k])/beta[k] # (10.61)
        W[,,k] = W0inv + Nk[k]*S[,,k] + ((beta0*Nk[k])/(beta0+Nk[k]))*tcrossprod(xbar[,k]-m0) # (10.62)
        W[,,k] = solve(W[,,k])
    }

    cat('xbar[,1]', xbar[,1], '\n')
    cat('S[,,1]', S[,,1], '\n')
    cat('m[,1]', m[,1], '\n')
    cat('W[,,k]', W[,,1], '\n')

    model$alpha = alpha
    model$m = m
    model$W = W
    model$v = v
    model$beta = beta
    model$S = S
    model$xbar = xbar
    model
}
