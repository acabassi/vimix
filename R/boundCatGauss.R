#' Lower bound for multiple categorical variables
#' @param X NxD matrix data
#' @param model Model parameters
#' @param prior Model parameters
#' @return Value of lower bound
#' @export
boundCatGauss = function(X, model, prior){

    alpha0 = prior$alpha
    eps0 = prior$eps

    alpha = model$alpha
    eps = model$eps # K x D x nCat

    Resp = model$Resp
    logResp = model$logResp

    N = dim(X)[1]
    D = dim(X)[2]
    K = dim(Resp)[2]
    nCat = dim(eps)[3]

    logpi = digamma(alpha) - digamma(sum(alpha))
    Epz = sum(Resp%*%logpi)
    # Resp[which(Resp<.Machine$double.xmin)] = .Machine$double.xmin # TODO controllare che questo non alteri i risultati
    Eqz = sum(Resp*logResp) # (10.75) univariate

    logCalpha0 = lgamma(K*alpha0)-K*lgamma(alpha0)
    Eppi = logCalpha0 + sum((alpha0-1)*logpi)# (10.73)

    logCalpha = lgamma(sum(alpha))-sum(lgamma(alpha))
    Eqpi = logCalpha + sum((alpha-1)*logpi)# (10.76)

    logCeps = array(0, c(K, D))
    ElnPhi = array(0, c(K,D,nCat))
    EpX = EpPhi = EqPhi = 0

    for(k in 1:K){ ## TODO this is extremely inefficient but it's just to check if everything works
        for(d in 1:D){

            maxXd = max(X[,d])
            logPhi0 = digamma(eps0[d,1:maxXd]) - digamma(sum(eps0[d,1:maxXd]))
            logCeps0 = lgamma(sum(eps0[d,1:maxXd])) - sum(lgamma(eps0[d,1:maxXd]))
            EpPhi = EpPhi + logCeps0 + sum((eps0[d,1:maxXd]-1)*logPhi0)

            logPhi = digamma(eps[k,d,1:maxXd]) - digamma(sum(eps[k,d,1:maxXd]))
            logCeps = lgamma(sum(eps[k,d,1:maxXd])) - sum(lgamma(eps[k,d,1:maxXd]))
            EqPhi = EqPhi + logCeps + sum((eps[k,d,1:maxXd]-1)*logPhi)

            for(n in 1:N){
                EpX = EpX + Resp[n,k]*logPhi[X[n,d]]
            }
        }
    }

    L = Epz - Eqz + Eppi - Eqpi + EpPhi - EqPhi + EpX

    if(!is.finite(L)) stop("Lower bound is not finite")

    L
}
