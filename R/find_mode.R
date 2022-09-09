find_mode <- function(fit){
  shape <- fit$shape_hat
  rate <- fit$rate_hat
  mode <- (shape - 1)/rate
  return(mode)
}