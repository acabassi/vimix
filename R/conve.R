#' Check convergence of variational algorithm
#' @param L Lower bound
#' @param iter Iteration number
#' @param tol Tolerance
#' @param maxiter Maximum number of iterations
#' @param verbose Boolean flag that, if TRUE, prints the outcome of the convergence check
#' @return Boolean flag that indicates if algorithm has converged
#' @export
check_convergence = function(L, iter, tol, maxiter, verbose){
    if (abs(L[iter]-L[iter-1]) < tol){ # stopping criterion
        if(verbose) message(sprintf("Converged in %d steps.\n", iter))
        conv = TRUE
    }else{
        conv = FALSE
        if(verbose) message(sprintf("Lower bound %f. \n", L[iter]))
    }
    if (iter == maxiter) warnings(sprintf("Not converged in %d steps.\n", maxiter))
    conv
}
