#' Expectation step for mixture of independent Gaussians
#'
#' @param X NxD data matrix.
#' @param model Model parameters.
#' @return Updated model parameters.
#' @export
#'
expectSelGauss = function(X, model, null){

    alpha = model$alpha
    v = model$v
    beta = model$beta
    m = model$m
    W = model$W
    c = model$c

    lnf_null = null$lnf

    N = dim(X)[1]
    D = dim(X)[2]
    K = length(alpha)

    logRho = matrix(NA, N, K)
    Elnf = array(NA, c(D, N, K))

    ElnPi = digamma(alpha) - digamma(sum(alpha))

    for (k in 1:K){

        ElnLa = D*log(2*pi) + sum(log(0.5 * (1/W[,k]))) + digamma(0.5*v[k]) # (10.65)
        diff = sweep(X, 2, m[,k], FUN="-")
        ExmuLaxmu = 1/beta[k]

        secondTerm = rep(0, N)
        for(d in 1:D){
            ExmuLaxmu = v[k]*(diff[,d]^2)*W[d,k] # (10.64)
            Elnf[d,,k] = 0.5*ElnLa - 0.5*ExmuLaxmu # (10.46)
            secondTerm = secondTerm + c[d]*Elnf[d,,k] + (1-c[d])*lnf_null[,d]
        }
        logRho[,k] = ElnPi[k] + secondTerm

    }

    logSumExpLogRho = apply(logRho, 1, log_sum_exp)

    logResp =  sweep(logRho, MARGIN = 1, STATS = logSumExpLogRho, FUN = "-")# 10.49
    Resp = apply(logResp, 2, exp)

    model$Elnf = Elnf
    model$logResp = logResp
    model$Resp = Resp
    model
}
