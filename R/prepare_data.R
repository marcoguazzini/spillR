#' Generate synthetic data according to the hierarchical model
#'
#' @export
#'
#' @param spillover_matrix Square matrix of spillover values
#'
#' @param n_cells Number of iterations
#' 
#' @param rates Gamma distribution rates
#' 
#' @param shapes Gamma distribution rates
#'
#' @param n_cells Number of iterations
#' 
#' @return matrix
#'
#' @examples
#' set.seed(23)
#' sm <- spillR::load_spillover()
#' shapes <- 10
#' rates <- 1
#' n_cells <- 5000
#' counts <- prepare_data(sm, shapes, rates, n_cells)


prepare_data <- function(spillover_matrix,shape,rate, n_cells){
 
    n <- ncol(spillover_matrix)
    channel_names <- rownames(spillover_matrix)
    lambdas = matrix(rgamma(n*n_cells, shape = shape, rate = rate),
                     nrow=n_cells, byrow=TRUE)
    counts = matrix(rpois(n*n_cells, lambda = lambdas %*% spillover_matrix),
                    nrow=n_cells, byrow=TRUE)
    colnames(counts) <- channel_names
    return(counts)
}
