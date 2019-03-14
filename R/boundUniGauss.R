#' Lower bound for univariate Gaussian
#' @param X NxD matrix data
#' @param model Model parameters
#' @param prior Model parameters
#' @return Value of lower bound
#' @export
boundUniGauss = function(X, model, prior){

    alpha0 = prior$alpha
    beta0 = prior$beta
    v0 = prior$v
    m0 = prior$m
    W0 = prior$W
    W0inv = prior$Winv

    alpha = model$alpha
    beta = model$beta
    v = model$v
    m = model$m
    W = model$W
    Resp = model$Resp
    logResp = model$logResp

    N = length(X)
    D = 1
    K = dim(Resp)[2]

    Nk = colSums(Resp)

    logpi = digamma(alpha) - digamma(sum(alpha))
    Epz = sum(Resp%*%logpi)
    # Resp[which(Resp<.Machine$double.xmin)] = .Machine$double.xmin # TO DO controllare che questo non alteri i risultati
    Eqz = sum(Resp*logResp) # (10.75) univariate

    logCalpha0 = lgamma(K*alpha0)-K*lgamma(alpha0)
    Eppi = logCalpha0 + sum((alpha0-1)*logpi)# (10.73)

    logCalpha = lgamma(sum(alpha))-sum(lgamma(alpha))
    Eqpi = logCalpha + sum((alpha-1)*logpi)# (10.76)

    EpX = EpMuLambda = EpMuLambda2 = EqMuLambda <- 0
    log_Lambda <- rep(0,K)

    for(k in 1:K){
        log_Lambda[k] <- digamma(v[k]/2) + log(0.5/W[k])
        EpX <- EpX + Nk[k] * (log_Lambda[k] - D/beta[k] -
                                  v[k]*W[k]*sum((X - m[k])^2) - D*log(2*pi))

        # (10.74)
        EpMuLambda <- EpMuLambda + D*log(beta0/(2*pi)) + log_Lambda[k] -
            (D*beta0)/beta[k] - beta0*v[k]*((m[k] - m0)^2) * W[k]
        EpMuLambda2 <- EpMuLambda2 + v[k] * W[k] / W0

        # (10.77)
        EqMuLambda <- EqMuLambda + 0.5*log_Lambda[k] + 0.5*D*log(beta[k]/(2*pi)) - logUniB(W[k], v[k]) -
            0.5*(v[k] - D - 1)*log_Lambda[k] + 0.5*v[k]*D
    }

    EpMuLambda <- 0.5*EpMuLambda + K*(logUniB(W0, v0)) + 0.5*(v0 - D - 1)*sum(log_Lambda) - 0.5*EpMuLambda2 # 10.74
    EpX  <- 0.5 * EpX # (10.71)

    L = Epz - Eqz + Eppi - Eqpi + EpMuLambda - EqMuLambda + EpX
    if(!is.finite(L)) stop("Lower bound is not finite")
    L
}

#' Compute univariate equivalent of logB
#' @param W Da rivedere
#' @param nu Da rivedere
logUniB <- function(W, nu){
    Winv <- 1/W
    if(is.na(lgamma(0.5 * nu))) stop()
    return(lgamma(0.5 * nu) + 0.5 * nu * log(0.5*Winv))
}
