#' Plot spillover curve and uncorrected/corrected counts
#'
#' @aliases plot.spillr
#' @method plot spillr
#'
#' @import ggplot2
#' @import magrittr
#' @import dplyr
#' @importFrom methods is
#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer
#' @importFrom cowplot plot_grid
#' @export
#'
#' @param x A \code{spillr} class
#' @param ... Other parameters
#' @return \code{\link[ggplot2]{ggplot2}} object
#'
#' @examples
#' set.seed(23)
#' tb_real <- generate_real()
#' tb_bead <- generate_bead()
#' target_marker <- "A"
#' spillover_marker <- "B"
#' comp <- spillR::compensate(tb_real, tb_bead, "A", "B")
#' plot(comp)
plot.spillr <- function(x, ...) {
  
  if(!is(x, "spillr"))
    stop("Input needs to be a spillr object computed by compensate function.")
  tfm <- function(x) asinh(x/5)
  tb_bead <- x$tb_bead
  tb_compensate <- x$tb_compensate
  tb_spill_prob <- x$tb_spill_prob
  target_marker <- x$target_marker
  
  p_spill <- 
    ggplot(tb_bead, aes(tfm(.data[[target_marker]]), color = barcode)) +
    geom_density(adjust = 1) +
    facet_wrap(~type, ncol = 1, scales = "free_y") +
    geom_line(data = tb_spill_prob, 
              aes(tfm(.data[[target_marker]]), spill_prob_smooth), 
              color = "black", 
              linetype = "longdash") + 
    ylab("density")
  
  p_comp <- tb_compensate %>% 
    pivot_longer(c(corrected, uncorrected), 
                 names_to = "status", values_to = target_marker) %>%
    ggplot(aes(tfm(.data[[target_marker]]))) + 
    geom_histogram(position = "identity", alpha = 0.5, bins = 50) + 
    facet_wrap(~status, nrow = 1) +
    geom_line(
      data = bind_rows(
        mutate(tb_spill_prob, status = "corrected"), 
        mutate(tb_spill_prob, status = "uncorrected")
      ) %>%
        mutate(spill_prob_smooth = spill_prob_smooth*500), 
      aes(tfm(.data[[target_marker]]), spill_prob_smooth), 
      color = "black", 
      linetype = "longdash") + 
    ylab("count")
  
  plot_grid(
    p_comp, 
    p_spill, 
    rel_widths = c(0.6, 0.4),
    nrow = 1)
  
}  
