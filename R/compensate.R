#' Generate compensated counts dataset
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import future.apply
#' @import extraDistr
#' @export
#'
#' @return \code{\link[tibble]{tibble}} data frame
#'
#' @examples
#' set.seed(23)
#' sm <- spillR::load_spillover()
#' counts <- spillR::prepare_data(sm, shape, rate, n_cells)
#' compensated_counts <- generate_data(counts,spillover_matrix)
#' compensated_counts





compensate <- function(counts,spillover_matrix){

plan(multisession, workers = 8)

  
  #lambda <- rep(0, length(channel_names))
  #names(lambda) <- channel_names
  #tb_info <- as_tibble(t(lambda))
  #tb_info <- tb_info[-1,]
  
  # lambda equal a zero vector, so that offset is zero
  channel_names <- colnames(counts)
  lambda <- rep(0, length(channel_names))
  names(lambda) <- channel_names
  # update fit including spillover by iterating over all markers
  n_iter <- 10 # this is the number of times to iterate over the graph
  modes <- matrix(NA, nrow = n_iter, ncol = length(channel_names))
  shapes <- matrix(NA, nrow = n_iter, ncol = length(channel_names))
  rates <- matrix(NA, nrow = n_iter, ncol = length(channel_names))
  colnames(modes) = colnames(shapes) = colnames(rates) <- channel_names
  for(iter in 1:n_iter) {
    
    update <- function(j) {
      # prepare and fit
      y <- counts[,j]                        # this is the response marker
      s <- spillover_matrix[-j,j]                         # s is incoming spillover
      offset <- sum(s*lambda[-j])            # sum of s_j lambda_j products
      fit <- fit_gampois(y, offset = offset) # same now with t and offset
      }
      # update current estimate of lambda, three options:
      # 1) mode
     # find_mode(fit)
      # 2) draw from gamma
      # rgamma(n = 1, shape = fit$shape_hat, rate = fit$rate_hat)
      # 3) MCMC style proposal
      # proposal <- rnorm(n = 1, mean = find_mode(fit), sd = 1)
      # if(proposal < 0)
      #   proposal <- 0
      # proposal

    #lambda_new <- future_sapply(seq(ncol(counts)), update, future.seed = TRUE)
    # update for next iteration
    #names(lambda_new) = channel_names
    #lambda <- lambda_new
    
    fit_list <- future_lapply(seq(ncol(counts)), update, future.seed = TRUE)
  
  # bookkeeping
  modes[iter,] <- sapply(fit_list, function(fit) find_mode(fit))
  shapes[iter,] <- sapply(fit_list, function(fit) fit$shape_hat)
  rates[iter,] <- sapply(fit_list, function(fit) fit$rate_hat)
  # update for next iteration
  lambda <- modes[iter,]
  }
                         
                         
                                            
   # generate with our method spillR::compensate
lambdas_spillr <- matrix(
  rgamma(n_cells*length(channel_names), 
         shape = shapes[n_iter,], 
         rate = rates[n_iter,]),
  nrow = n_cells,
  ncol = length(channel_names), 
  byrow = TRUE
  )
compensated_counts <- matrix(
  rpois(length(lambdas_spillr), lambda = lambdas_spillr),
  nrow = n_cells,
  ncol = length(channel_names)
  )
  colnames(compensated_counts) <- channel_names
 
  return(compensated_counts)
}
