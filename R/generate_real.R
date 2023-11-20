#' Generate dataset for vignettes and simulation studies
#'
#' @importFrom tibble tibble
#' @export
#'
#' @return \code{\link[tibble]{tibble}} data frame
#'
#' @examples
#' set.seed(23)
#' generate_real()
generate_real <- function(){
    n_real <- 10000
    lambda_real <- 100
    lambda_bead  <- 70
    spill_prob   <- 0.1
    # real experiment
    Z_target <- rpois(n=n_real, lambda=lambda_real)
    Z_spill  <- rpois(n=n_real, lambda=lambda_bead)
    spill    <- rbinom(n=n_real, size=1, prob=spill_prob)
    Y        <- (1-spill)*Z_target + spill*Z_spill
    df_real  <- tibble(Y=Y, barcode="Y", type="real cells")
    df_real
    
}
