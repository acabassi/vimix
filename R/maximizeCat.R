#' Maximization step for categorical variables
#'
#' @param X NxD data matrix
#' @param model Model parameters
#' @param prior Prior parameters
#' @return Updated model parameters
#' @export
#'
maximizeCatGauss = function(X, model, prior){

    alpha0 = prior$alpha
    eps0 = prior$eps
    Resp = model$Resp

    eps = model$eps

    N = dim(X)[1]
    D = dim(X)[2]
    K = dim(Resp)[2]
    maxNCat = dim(eps)[3]

    Nk = colSums(Resp) + 1e-10 # (10.51)
    alpha = alpha0 + Nk # (10.58)

    eps = array(0, c(K, D, maxNCat))
    for(d in 1:D){
        for(j in 1:maxNCat){
            eps[,d,j] = eps0[d,j] + t(Resp)%*%(X[,d]==j)
        }
    }

    model$alpha = alpha
    model$eps = eps
    model
}
