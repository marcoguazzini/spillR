prepare_data <- function(spillover_matrix,shape,rate, n_cells){
    #' Generate synthetic data according to the hierarchical model
    #' Args:
    #'     spillover_matrix: square matrix of spillover values
    #'     n_sample: number of iterations
    #'     rates: gamma distribution rates
    #'     shapes: gamma distribution shapes
    n <- ncol(spillover_matrix)
    channel_names <- rownames(spillover_matrix)
    lambdas = matrix(rgamma(n*n_cells, shape = shape, rate = rate),
                     nrow=n_cells, byrow=TRUE)
    counts = matrix(rpois(n*n_cells, lambda = lambdas %*% spillover_matrix),
                    nrow=n_cells, byrow=TRUE)
    colnames(counts) <- channel_names
    return(counts)
}
