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
                        starting_value = c(1, 0.1, 0)) {
  dgpoisshifted_zeroinf <- function(y, r, b, p, t, log = FALSE) {
    
    if (!log) {
      factor1 <- exp(-t) * b ^ r / (factorial(y) * gamma(r))
      term <-
        function(k)
          choose(y, k) * t ^ (y - k) * gamma(k + r) / (b + 1) ^ (k + r)
      factor2 <- sum(sapply(0:y, term))
      return(p * (y == 0) + (1 - p) * (factor1 * factor2))
    } else {
      if (y > 0) {
        factor1 <- -t + r * log(b) - (lfactorial(y) + lgamma(r))
        term <-
          function(k)
            lchoose(y, k) + (y - k) * log(t) + lgamma(k + r) - (k + r) * log(b + 1)
        factor2 <- 0
        if (t > 0) {
          factor2 <- logSumExp(sapply(0:y, term))
        } else {
          factor2 <- lgamma(y + r) - (y + r) * log(b + 1)
        }
        # will produce -Inf if p = 1, 
        # but p = 1 means that there are only zeros,
        # we exclude 1 by setting the limit of p to 1 - epsilon
        return(log(1 - p) + factor1 + factor2)
      }
      else {
        lb <- log(1 - p) - t + r * log(b) - r * log(b + 1)
        if(p == 0) {
          lx <- lb
        } else {
          la <- log(p)
          lx <- logSumExp(c(la, lb))
        }
        return(lx)
      }
    }
  }
  
  likelihood <- function(param, x) {
    shape <- param[1]
    rate <- param[2]
    p <- param[3]
    ll <- dgpoisshifted_zeroinf(x, shape, rate, p, t = offset, log = TRUE)
    return(ll)
  }
  log_lik <- function(param, x) {
    sumll <- -sum(sapply(x, function(xi)
      likelihood(param, xi)))
    return(sumll)
  }
  res <- optim(
    starting_value,
    log_lik,
    x = x,
    method = "L-BFGS-B",
    lower = c(.Machine$double.eps, .Machine$double.eps, 0),
    upper = c(Inf, Inf, 1-.Machine$double.eps)
  )$par
  
  data.frame(shape_hat = res[1],
             rate_hat = res[2],
             p_hat = res[3])
}
