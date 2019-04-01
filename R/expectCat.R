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
