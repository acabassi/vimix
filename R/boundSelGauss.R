#' Lower bound for independent Gaussian
#' @param X NxD matrix data
#' @param model Model parameters
#' @param prior Model parameters
#' @return Value of lower bound
#' @export
#'
boundSelGauss = function(X, model, prior, null){

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
    S = model$S
    xbar = model$xbar
    Resp = model$Resp
    logResp = model$logResp

    m_null = null$m_null
    sd_null = null$stdev

    N = dim(X)[1]
    D = dim(X)[2]
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

    EpX = EpMuLambda1 = EpMuLambda2 = EqMuLambda <- 0
    log_Lambda <- rep(0,K)

    for(k in 1:K){
        log_Lambda[k] <- D*digamma(v[k]/2) + sum(log(0.5 /W[,k]))
        pollo <- sweep(X, 2, m[,k], FUN = "-")
        EpX <- EpX + Nk[k] * (log_Lambda[k] - D/beta[k] -
                                  v[k]*sum((pollo)^2%*%W[,k]) - D*log(2*pi))

        # (10.74)
        EpMuLambda1 <- EpMuLambda1 + D*log(beta0/(2*pi)) + log_Lambda[k] -
            (D*beta0)/beta[k] - sum(beta0*v[k]*((m[,k] - m0)^2) * W[,k])
        EpMuLambda2 <- EpMuLambda2 + sum(v[k] * W[,k] / W0)

        # (10.77)
        EqMuLambda <- EqMuLambda + 0.5*D*log(beta[k]/(2*pi)) - sum(logUniB(W[,k], v[k])) -
            D*0.5*(v[k] - 1)*log_Lambda[k] + 0.5*D*(v[k]+1)
        # I'm not 100% sure about the D*(v0 *0.5 - 1)*sum(log_Lambda) term.
        # It would be nicer if it matched the multivariate case: 0.5*(v[k] - D - 1)*log_Lambda[k]
    }

    EpMuLambda <- 0.5*EpMuLambda1 + K*D*sum(logUniB(W0, v0)) + D*(v0 *0.5 - 1)*sum(log_Lambda) - 0.5*EpMuLambda2 # 10.74
    # I'm not 100% sure about the D*(v0 *0.5 - 1)*sum(log_Lambda) term.
    # It would be nicer if it matched the multivariate case:  0.5*(v0 - D - 1)*sum(log_Lambda)

    EpX  <- 0.5 * EpX # (10.71)

    L = Epz - Eqz + Eppi - Eqpi + EpMuLambda - EqMuLambda + EpX

    if(!is.finite(L)) stop("Lower bound is not finite")

    L
}
