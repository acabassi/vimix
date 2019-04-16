#' Expectation step for mixture of independent Gaussians
#'
#' @param X NxD data matrix.
#' @param model Model parameters.
#' @param prior Prior hyperparameters.
#' @return Updated model parameters.
#' @export
#'
updateSalGauss = function(X, model, prior){

    alpha = model$alpha
    beta = model$beta
    m = model$m
    c = model$c
    pi = model$pi
    w = model$w
    gamma = model$gamma
    Resp = model$Resp
    eps = model$eps
    xi = model$xi
    rho = model$rho

    c0 = prior$c
    alpha0 = prior$alpha
    beta0 = prior$beta
    m0 = prior$m
    
    N = dim(X)[1]
    D = dim(X)[2]
    K = length(pi)
 
    rtilde = matrix(NA, N, K)
    
    for(n in 1:N){
        for(k in 1:K){
            temp = 0
            for(d in 1:D){
                temp =  temp + 0.5* rho[n,d] * ( digamma(alpha[d,k]) - log(beta[d,k]) - 
                                                     alpha[d,k]*((X[n,d]-m[d,k])^2 + 1/c[d,k])/beta[d,k] )
            }
            rtilde[n,k] = exp(temp)
        }
    }
    
    for(n in 1:N){
        for(k in 1:K){
            Resp[n,k] = pi[k] * rtilde[n,k] / sum( pi * rtilde[n,] )
            
        }
    }
    
    for(d in 1:D){
        for(k in 1:K){
            m_numer = c0 * m0[d] + (alpha[d,k] / beta[d,k]) * sum(Resp[,k] * rho[,d] * X[,d])
            m_denom = c0 + (alpha[d,k] / beta[d,k]) * sum(Resp[,k] * rho[,d])
            m[d,k] = m_numer / m_denom
            
            c[d,k] = c0 + alpha[d,k] * sum( Resp[,k] * rho[,d] ) / beta[d,k]
            
            alpha[d,k] = alpha0 + 0.5 * sum( Resp[,k] * rho[,d] )
            
            beta[d,k] = beta0 + 0.5 * sum (Resp[,k] * rho[,d] * ((X[,d] - m[d,k])^2 + 1/c[d,k]) )
            
        }
    }
    
    for(d in 1:d){
        for(n in 1:N){
            temp = 0
            for(k in 1:K){
                temp = temp + 0.5 * Resp[n,k] * ( digamma(alpha[d,k]) - log(beta[d,k]) - 
                                                      alpha[d,k] * ((X[n,d] -  m[d,k])^2 + 1/c[d,k])  / beta[d,k] )
            }
            rhotilde[n,d] = exp(temp)
            xi[n,d] = exp(0.5 * (- gamma[d] * ( X[n,d] - eps[d] )^2 + log( gamma[d] )))
            rho[n,d] = w[d] * rhotilde[n,d] / (w[d] * rhotilde[n,d])
            
        }
        
        w[d] = sum( rho[,d] ) / N
        eps[d] = sum( rho[,d] * X[,d] ) / sum(rho[,d])
        gamma[d] = sum( rho[,d] ) / sum( rho[,d] * ( X[,d] - eps[d] )^2 )
        
    }
    
    for(k in 1:K){
        pi[k] = sum(Resp[,k])/N
    }
 
    model = list(alpha = alpha, beta = beta, m = m, c = c, pi = pi, w = w, 
                 gamma = gamma, Resp = Resp + 1e-10, eps = eps, xi = xi)
}
