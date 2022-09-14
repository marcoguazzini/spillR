find_mode <- function(fit) {
  mode <- 0
  if(fit$shape_hat >= 1) {
    mode <- (fit$shape_hat-1)/fit$rate_hat
  }
  mode
}
