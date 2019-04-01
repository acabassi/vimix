initializeC = function(X, model, prior, null){

    alpha = model$alpha
    v = model$v
    beta = model$beta
    m = model$m
    W = model$W

    m_null = null$m_null
    W_null = null$W_null

    N = dim(X)[1]
    D = dim(X)[2]
    K = length(v)

    C = rep(0, D)

    for(d in 1:D){

        # This may not be necessary
        # ElnGammad = digamma()
        # nu1d =
    }

}
