compensate <- function(counts,spillover_matrix){
library(tibble)
library(dplyr)
library(tidyr)
library(extraDistr)
library(future.apply)
plan(multisession, workers = 8)
  # lambda equal a zero vector, so that offset is zero
  smc_scaled <- 10*spillover_matrix
  diag(smc_scaled) <- 1.0
  channel_names <- colnames(counts)
  lambda <- rep(0, length(channel_names))
  names(lambda) <- channel_names
  tb_info <- as_tibble(t(lambda))
  tb_info <- tb_info[-1,]
  # update fit including spillover by iterating over all markers
  n_iter <- 10 # this is the number of times to iterate over the graph
  for(iter in 1:n_iter) {
    
    update <- function(j) {
      # prepare and fit
      y <- counts[,j]                        # this is the response marker
      s <- smc_scaled[-j,j]                         # s is incoming spillover
      offset <- sum(s*lambda[-j])            # sum of s_j lambda_j products
      fit <- fit_gampois(y, offset = offset) # same now with t and offset
      
      # update current estimate of lambda, three options:
      # 1) mode
      find_mode(fit)
      # 2) draw from gamma
      # rgamma(n = 1, shape = fit$shape_hat, rate = fit$rate_hat)
      # 3) MCMC style proposal
      # proposal <- rnorm(n = 1, mean = find_mode(fit), sd = 1)
      # if(proposal < 0)
      #   proposal <- 0
      # proposal
    }
    lambda_new <- future_sapply(seq(ncol(counts)), update, future.seed = TRUE)
    # update for next iteration
    names(lambda_new) = channel_names
    lambda <- lambda_new
    
    # bookkeeping
    tb_info <- bind_rows(tb_info, lambda)
  }
  tb_info <- mutate(tb_info, n = 1:nrow(tb_info))
  n = ncol(spillover_matrix)
  n_sample <- nrow(counts)
  val <- as.numeric(tb_info[nrow(tb_info),1:n])
  lambda_compensated = matrix(val, n_cells, ncol=ncol(tb_info)-1, byrow=TRUE)
  compensated_counts = matrix(rpois(n*n_sample, lambda = lambda_compensated),
                              nrow = n_sample, byrow=TRUE)
  colnames(compensated_counts) <- channel_names
 
  return(compensated_counts)
}
