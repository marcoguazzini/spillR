#' Compute the mode of a gamma distibuted random variabl
#'

#' @param fit Vector with the two gamma's parameter shape and rate
#' 
#' @return value
#'
#' @examples
#' set.seed(23)
#' shape <- 10
#' rate <- 1
#' param <- c(shape,rate)
#' find_mode(param)



find_mode <- function(fit) {
  mode <- 0
  if(fit$shape_hat >= 1) {
    mode <- (fit$shape_hat-1)/fit$rate_hat
  }
  mode
}
