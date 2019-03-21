#' Preprocessing matrix containing categorical variables
#'
#' @param X NxD matrix
#' @references Inspiration came from code contained in the R package 'mixdir'
#' @export
#'
preprocessCategoricalDataset = function(X){
    # Create a numeric matrix with the entries 1 to N_cat_j
    X <- as.data.frame(X)

    X[colnames(X)] <- lapply(X[colnames(X)], function(col){
        if(! is.factor(col)){
            col <- as.factor(as.character(col))
        }
        levels(col) <- c(levels(col), "(Missing)")
        col[is.na(col)] <- "(Missing)"
        col
    })
    categories <- lapply(1:ncol(X), function(j)levels(X[, j]))

    nCat <- sapply(categories, length) # make sure that there are no empty categories
    if(any(nCat == 0)) {
        stop("Column ", paste0(colnames(X)[nCat == 0], collapse = ","), " is empty (i.e. all values are NA). Please fix this.")
    }

    X[colnames(X)] <- lapply(X[colnames(X)], function(col) as.numeric(col))
    X <- as.matrix(X)
}
