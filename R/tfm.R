#' Variance stabilizing transform of counts
#'
#' @param x Raw count
#'
#' @return A transformed count
tfm <- function(x) {
    asinh(x / 5)
}
