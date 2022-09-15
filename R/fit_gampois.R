#' Perform a Maximum Likelihood Estimation of a shifted gamma Poisson distribution
#'
#' @importFrom extraDistr dgpois
#' @import matrixStats
#' @export
#'
#' @param x Vector of values representing the samples
#'
#' @param offset Value to be added in order to shift the gamma distribution
#'  
#' @param starting_value Starting value for the optimization matrix
#'
#' @return data_frame
#'
#' @examples
#' set.seed(23)
#' n_sample <- 1000
#' shape <- 10
#' rate <- 1
#' lambda <- rgamma(n_sample, shape, rate)
#' offset <- 10
#' x <- rpois(n_sample, lambda = lambda + offset)
#' fit <- fit_gampois(x = x, offset = offset, starting_value = c(1,0.1))
#' fit





fit_gampois <- function(x,
                        offset = 0,
                        starting_value = c(1,0.1)) {
  dgpoisshifted <- function(y, r, b, t, log = FALSE) {
    pm <- 0
    if(!log) {
      factor1 <- exp(-t)*b^r / (factorial(y)*gamma(r))
      term <- function(k) choose(y, k)*t^(y-k)*gamma(k+r)/(b+1)^(k+r)
      factor2 <- sum(sapply(0:y, term))
      pm <- factor1 * factor2
    } else {
      factor1 <- -t+log(b)*r-(lfactorial(y)+lgamma(r))
      term <- function(k) lchoose(y,k)+(y-k)*log(t)+ lgamma(k+r)-(k+r)*log(b+1)
      factor2 <- logSumExp(sapply(0:y, term))
      pm <- factor1 + factor2
    }
    pm
  }
  likelihood <- function(param,x){
    shape <- param[1]
    rate <- param[2]
    if(offset == 0){
      ll <- dgpois(x = x, 
                   shape = shape , 
                   rate = rate,
                   log = TRUE)
    }
    else{
      ll <- dgpoisshifted(x,shape,rate, t = offset,log=TRUE)}
    return(ll)
  }
  log_lik <- function(param, x) {
    sumll <- -sum(sapply(x, function(xi)
      likelihood(param, xi)))
    return(sumll)
  }
  res <- optim(
    starting_value,
    log_lik ,
    x = x,
    method = "L-BFGS-B",
    lower = c(1, 0.001)
  )$par
  
  data.frame(shape_hat = res[1], 
             rate_hat = res[2]) 
}
