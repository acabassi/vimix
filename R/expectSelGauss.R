#' Expectation step for mixture of independent Gaussians
#'
#' @param X NxD data matrix.
#' @param model Model parameters.
#' @return Updated model parameters.
#' @export
#'
expectSelGauss = function(X, model){

    alpha = model$alpha
    v = model$v
    beta = model$beta
    m = model$m
    W = model$W
    cc = model$c

    lnf_null = model$lnf

    N = dim(X)[1]
    D = dim(X)[2]
    K = length(alpha)

    logRho = matrix(NA, N, K)
    Elnf = array(NA, c(D, N, K))

    ElnPi = digamma(alpha) - digamma(sum(alpha))

    for (k in 1:K){

        diff = sweep(X, 2, m[,k], FUN="-")

        secondTerm = rep(0, N)
        for(d in 1:D){
            Elnf[d,,k] = 0.5 * (log(2*pi) - log(0.5/W[d,k]) + digamma(0.5*v[d,k]) - v[d,k]*(diff[,d]^2)*W[d,k] - 1/beta[d,k]) # (10.46)
            secondTerm = secondTerm + cc[d]*Elnf[d,,k] + (1-cc[d])*lnf_null[,d]
        }
        logRho[,k] = ElnPi[k] + secondTerm

    }

    logSumExpLogRho = apply(logRho, 1, log_sum_exp)

    logResp =  sweep(logRho, MARGIN = 1, STATS = logSumExpLogRho, FUN = "-")# 10.49
    Resp = apply(logResp, 2, exp)

    model$Elnf = Elnf
    model$logResp = logResp
    model$Resp = Resp + 1e-10
    model
}
