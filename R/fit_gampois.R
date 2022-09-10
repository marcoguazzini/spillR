fit_gampois <- function(x,
                        offset = 0,
                        starting_value = c(1,0.1)) {
  library(extraDistr)
  library(matrixStats)
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
