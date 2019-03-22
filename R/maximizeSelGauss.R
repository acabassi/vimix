#' Maximization step for categorical variables
#'
#' @param X NxD data matrix
#' @param model Model parameters
#' @param prior Prior parameters
#' @return Updated model parameters
#' @export
#'
maximizeSelGauss = function(X, model, prior, null){

    alpha0 = prior$alpha
    beta0 = prior$beta
    m0 = prior$m
    v0 = prior$v
    W0 = prior$W
    o = prior$o

    Resp = model$Resp
    c = model$c
    Elnf = model$Elnf # DxNxK

    lnf_null = null$lnf

    N = dim(X)[1]
    D = dim(X)[2]
    K = dim(Resp)[2]

    xbar = m = W = S = matrix(0, D, K)

    c_new = rep(NA, D)
    for(d in 1:D){
        ElnGammad = digamma(c[d] + o) - digamma(2*d + 1)
        lnNu1d = ElnGammad + sum(Resp*Elnf[d,,])

        ElnOneMinusGammad = digamma(1-c[d]+o) - digamma(2*d + 1)
        lnNu2d = ElnOneMinusGammad + sum(Resp*matrix(lnf_null[,d], N, K, byrow = F))

        c_new[d] = lnNu1d / (lnNu1d + lnNu2d)
    }


    Nk = colSums(Resp) + 1e-10 # (10.51)
    Ndk = matrix(rep(Nk,D), K, D) # KxD matrix
    Ndk = t(Ndk * matrix(c, K, D, byrow = T))

    alpha = alpha0 + Nk # (10.58)
    beta = beta0 + Ndk # (10.60)
    v = v0 + Ndk # (10.63)

    for(k in 1:K){
        xbar[,k] = (Resp[,k]%*%X)/Nk[k] # (10.52)
        x_cen = sweep(X, MARGIN = 2, STATS = xbar[,k], FUN = "-")
        S[,k] = (t(x_cen^2)%*%Resp[,k])/Nk[k] # (10.53)
        for(d in 1:D){
            m[d,k] = (beta0*m0+Ndk[d,k]*xbar[d,k])/beta[k] # (10.61)
            W[d,k] = 1/W0[d] + Ndk[d,k]*S[d,k] + ((beta0*Ndk[d,k])/(beta0+Ndk[d,k]))*((xbar[d,k]-m0)^2) # (10.62)
            W[d,k] = 1/W[d,k]
        }
    }

    model$alpha = alpha
    model$m = m
    model$W = W
    model$v = v
    model$beta = beta
    model$S = S
    model$xbar = xbar
    model$c = c_new
    model
}
