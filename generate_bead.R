#' Generate dataset for vignettes and simulation studies
#'
#' @import tibble
#' @export
#'
#' @return \code{\link[tibble]{tibble}} data frame
#'
#' @examples
#' set.seed(23)
#' df <- generate_bead()
#' df
generate_bead <- function(){
  n_bead <- 1000
  lambda_bead  <- 70
  # real experiment
  Z_bead   <- rpois(n = n_bead, lambda = lambda_bead)
  df_bead  <- tibble(Y = Z_bead, barcode = "Z", type = "beads")
  df_bead
}
