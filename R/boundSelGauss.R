#' Lower bound for independent Gaussian
#' @param X NxD matrix data
#' @param model Model parameters
#' @param prior Model parameters
#' @return Value of lower bound
#' @export
#'
boundSelGauss = function(X, model, prior){

    alpha0 = prior$alpha
    beta0 = prior$beta
    v0 = prior$v
    m0 = prior$m
    W0 = prior$W
    W0inv = prior$Winv
    o = prior$o

    alpha = model$alpha
    beta = model$beta
    v = model$v
    m = model$m
    W = model$W
    S = model$S
    xbar = model$xbar
    Resp = model$Resp
    logResp = model$logResp
    Nk = model$Nk
    c = model$c

    lnf_null = model$lnf
    
    N = dim(X)[1]
    D = dim(X)[2]
    K = dim(Resp)[2]

    logpi = digamma(alpha) - digamma(sum(alpha))
    Epz = sum(Resp%*%logpi)
    # Resp[which(Resp<.Machine$double.xmin)] = .Machine$double.xmin # TO DO controllare che questo non alteri i risultati
    Eqz = sum(Resp*logResp) # (10.75) univariate

    logCalpha0 = lgamma(K*alpha0)-K*lgamma(alpha0)
    Eppi = logCalpha0 + sum((alpha0-1)*logpi)# (10.73)

    logCalpha = lgamma(sum(alpha))-sum(lgamma(alpha))
    Eqpi = logCalpha + sum((alpha-1)*logpi)# (10.76)

    EpX = EpMuLambda1 = EpMuLambda2 = EqMuLambda = EpGamma = EqGamma = EpDelta = EqDelta<- 0
    log_Lambda <- matrix(0, D, K)

    Elnf = array(NA, c(D, N, K))

    for(k in 1:K){

        pollo <- sweep(X, 2, m[,k], FUN = "-")

        for(d in 1:D){

            log_Lambda[d,k] <- digamma(0.5 * v[d,k]) - log(0.5 / W[d,k])

            # (10.74)
            EpMuLambda1 <- EpMuLambda1 + log(beta0/(2*pi)) -
                beta0/beta[d,k] - beta0*v[d,k]*((m[d,k] - m0[d])^2) * W[d,k]
            EpMuLambda2 <- EpMuLambda2 + v[d,k] * W[d,k] / W0[d]

            # (10.77)
            EqMuLambda <- EqMuLambda + 0.5*log(beta[d,k]/(2*pi)) - logUniB(W[d,k], v[d,k]) +
                0.5*(v[d,k] - 1)*log_Lambda[d,k] - 0.5*(v[d,k]+1)
        
            Elnf[d,,k] <- log_Lambda[d,k] - 1/beta[d,k] - v[d,k]*W[d,k]*sum((pollo[,d])^2) - log(2*pi)
            EpX <- EpX + ((c[d] * sum(Elnf[d,,k] * Resp[,k])) + (1-c[d]) * sum(lnf_null[,d] * Resp[,k]))

        }
    }
    
    EpMuLambda <- 0.5 * EpMuLambda1 + K * D * sum(logUniB(W0, v0)) + 0.5 * (v0 - 1) * sum(log_Lambda) - 0.5 * EpMuLambda2 # 10.74
    
    for(i in 1:D){
        ElnGammad <- digamma(c[d] + o) - digamma(2*o + 1)
        ElnOneMinusGammad <- digamma(1 - c[d] + o) - digamma(2*o + 1)
        EpGamma <- EpGamma + c[d] * ElnGammad + (1 - c[d]) * ElnOneMinusGammad
        
        lnNu1d = ElnGammad + sum(Resp*Elnf[d,,])
        lnNu2d = ElnOneMinusGammad + sum(Resp*matrix(lnf_null[,d], N, K, byrow = F))
        EqGamma <- c[d] * lnNu1d + (1 - c[d]) * lnNu2d
        
        EpDelta <- EpDelta + (o-1)* ElnGammad + (o-1) * ElnOneMinusGammad
        EqDelta <- EqDelta - lbeta(c[d]+o, 1-c[d]+o) + (c[d]+o-1)*ElnGammad + (o-c[d])*ElnOneMinusGammad
    }

    EpDelta <- EpDelta - D*lbeta(o,o)

    L = Epz - Eqz + Eppi - Eqpi + EpGamma - EqGamma + EpDelta - EqDelta + EpMuLambda - EqMuLambda + EpX

    if(!is.finite(L)) stop("Lower bound is not finite")

    return(list(L=L, Epz=Epz, Eqz=Eqz, Eppi=Eppi, Eqpi=Eqpi, EpGamma=EpGamma, EqGamma=EqGamma, 
           EpDelta=EpDelta, EqDelta=EqDelta, EpMuLambda=EpMuLambda, EqMuLambda = EqMuLambda, EpX=EpX))
    # L
}
