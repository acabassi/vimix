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

        ElnLa = D*log(2*pi) - log(0.5 * (1/W[k])) + digamma(0.5*v[k]) # (10.65)
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
