#' Expectation step for mixture of categorical variables
#'
#' @param X NxD data matrix.
#' @param model Model parameters.
#' @return Updated model parameters.
#' @export
#'
expectCatGauss = function(X, model){

    alpha = model$alpha
    eps = model$eps # K x D x nCat

    N = dim(X)[1]
    D = dim(X)[2]
    K = length(alpha)
    nCat = dim(eps)[3]

    ElnPi = digamma(alpha) - digamma(sum(alpha)) # K-dimensional vector
    logRho = matrix(ElnPi, N, K, byrow = TRUE)

    for(obs in 1:N){ ## TODO make this more efficient by using apply()
        for(k in 1:K){
            for(d in 1:D){
                maxXd = max(X[,d])
                logRho[obs,k] = logRho[obs,k] + digamma(eps[k,d,X[obs,d]]) - digamma(sum(eps[k,d,1:maxXd]))
            }
        }
    }

    logSumExpLogRho = apply(logRho, 1, log_sum_exp)

    logResp =  sweep(logRho, MARGIN = 1, STATS = logSumExpLogRho, FUN = "-")# 10.49
    Resp = apply(logResp, 2, exp)

    model$logResp = logResp
    model$Resp = Resp
    model
}

#' Expectation step for mixture of independent Gaussians
#'
#' @param X NxD data matrix.
#' @param model Model parameters.
#' @return Updated model parameters.
#' @export
#'
expectIndGauss = function(X, model){

    alpha = model$alpha
    v = model$v
    beta = model$beta
    m = model$m
    W = model$W

    N = dim(X)[1]
    D = dim(X)[2]
    K = length(v)

    logRho = matrix(NA, N, K)

    ElnPi = digamma(alpha) - digamma(sum(alpha))

    for (k in 1:K){

        ElnLa = D*log(2*pi) + sum(log(0.5 * (1/W[,k]))) + digamma(0.5*v[k]) # (10.65)
        diff = sweep(X, 2, m[,k], FUN="-")
        ExmuLaxmu = 1/beta[k]; for(d in 1:D) ExmuLaxmu = ExmuLaxmu + v[k]*(diff[,d]^2)*W[d,k] # (10.64)
        logRho[,k] = ElnPi[k] + 0.5*ElnLa - 0.5*ExmuLaxmu # (10.46)

    }

    logSumExpLogRho = apply(logRho, 1, log_sum_exp)

    logResp =  sweep(logRho, MARGIN = 1, STATS = logSumExpLogRho, FUN = "-")# 10.49
    Resp = apply(logResp, 2, exp)

    model$logResp = logResp
    model$Resp = Resp
    model
}

#' Expectation step for mixture of univariate Gaussians
#'
#' @param X NxD data matrix.
#' @param model Model parameters.
#' @return Updated model parameters.
#' @export
#'
expectUniGauss = function(X, model){

    alpha = model$alpha
    v = model$v
    beta = model$beta
    m = model$m
    W = model$W

    N = length(X)
    D = 1
    K = length(v)

    logRho = matrix(NA, N, K)

    ElnPi = digamma(alpha) - digamma(sum(alpha))

    for (k in 1:K){

        ElnLa = D*log(2*pi) + log(0.5 * (1/W[k])) + digamma(0.5*v[k]) # (10.65)
        ExmuLaxmu = 1/beta[k] + v[k]*((X-m[k])^2)*W[k] # (10.64)
        logRho[,k] = ElnPi[k] + 0.5*ElnLa - 0.5*ExmuLaxmu # (10.46)

    }

    logSumExpLogRho = apply(logRho, 1, log_sum_exp)

    logResp =  sweep(logRho, MARGIN = 1, STATS = logSumExpLogRho, FUN = "-")# 10.49
    Resp = apply(logResp, 2, exp)

    model$logResp = logResp
    model$Resp = Resp
    model
}

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
    W = model$W

    N = dim(X)[1]
    D = dim(X)[2]
    K = length(v)

    logRho = matrix(NA, N, K)

    ElnPi = digamma(alpha) - digamma(sum(alpha))

    for (k in 1:K){

        ElnLa = D*log(2) + log(det(W[,,k])) + sum(digamma(0.5*(v[k]+1-1:D))) # (10.65)

        diff = sweep(X, 2, m[,k], FUN="-")
        xmuLaxmu = diag(diff%*%W[,,k]%*%t(diff))
        ExmuLaxmu = D/beta[k] + v[k]*xmuLaxmu # (10.64)

        logRho[,k] = ElnPi[k] + 0.5*ElnLa - 0.5*ExmuLaxmu # (10.46)

    }


    logSumExpLogRho = apply(logRho, 1, log_sum_exp)

    logResp =  sweep(logRho, MARGIN = 1, STATS = logSumExpLogRho, FUN = "-")# 10.49
    Resp = apply(logResp, 2, exp)

    model$logResp = logResp
    model$Resp = Resp
    model
}

#' Use the log sum exp trick to improve numerical stability
#' @param x Vector of values
log_sum_exp <- function(x) {
    # Computes log(sum(exp(x))
    offset <- max(x)
    s <- log(sum(exp(x - offset))) + offset
    i <- which(!is.finite(s))
    if (length(i) > 0) { s[i] <- offset }
    return(s)
}
