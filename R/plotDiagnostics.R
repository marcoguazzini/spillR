#' Compute spillover probability and correct for spillover
#'
#' @importFrom CATALYST prepData assignPrelim applyCutoffs computeSpillmat
#'                      prepData
#' @importFrom dplyr as_tibble pull mutate
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot geom_density geom_freqpoly scale_color_manual
#'                     scale_linetype_manual geom_line xlab ylab ggtitle aes
#' @importFrom magrittr %>% %<>%
#' @importFrom SummarizedExperiment assay rowData assay<-
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#' @param ch Character string specifying the channel to plot
#'
#' @return A list of \code{\link[ggplot2]{ggplot2}} plots
#'
#' @examples
#' library(dplyr)
#' library(magrittr)
#' bc_key <- c(139, 141:156, 158:176)
#' sce_bead <- CATALYST::prepData(ss_exp)
#' sce_bead <- CATALYST::assignPrelim(sce_bead, bc_key, verbose = FALSE)
#' sce_bead <- CATALYST::applyCutoffs(estCutoffs(sce_bead))
#' sce_bead <- CATALYST::computeSpillmat(sce_bead)
#' data(mp_cells, package = "CATALYST")
#' sce <- CATALYST::prepData(mp_cells)
#' marker_to_barc <- rowData(sce_bead)[, c("channel_name", "is_bc")] %>%
#'   as_tibble() %>%
#'   filter(is_bc == TRUE) %>%
#'   mutate(barcode = bc_key) %>%
#'   select(marker = channel_name, barcode)
#' spillR::compCytof(sce, sce_bead, marker_to_barc, overwrite = FALSE)
#' plotDiagnostics(sce, "Yb173Di")
plotDiagnostics <- function(sce, ch) {
  tfm <- function(x) asinh(x / 5)

  # before and after correction
  exprs <- sce %>%
    assay("exprs") %>%
    t() %>%
    as_tibble()
  names(exprs) <- rowData(sce)$channel_name

  compexprs <- sce %>%
    assay("compexprs") %>%
    t() %>%
    as_tibble()
  names(compexprs) <- rowData(sce)$channel_name

  before_after <- tibble(
    cell = seq_along(nrow(exprs)),
    before = exprs %>% pull(ch),
    after = compexprs %>% pull(ch)
  )
  p_before_after <- before_after %>%
    pivot_longer(-.data$cell, names_to = "correction") %>%
    mutate(correction = factor(correction, levels = c("after", "before"))) %>%
    ggplot(aes(value, color = correction, linetype = correction)) +
    geom_freqpoly(alpha = 1.0, bins = 50, linewidth = 0.8) +
    scale_color_manual(values = c("#00BFC4", "#F8766D")) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    xlab(paste0("tfm(", ch, ")")) +
    ggtitle("Spillover Compensation on Real Cells")

  # diagnostic plot for spillover estimate
  tb_bead <- metadata(sce)$beads_distr[[ch]]
  tb_bead <- mutate(
    tb_bead,
    barcode = ifelse(barcode == ch, paste(ch, "(target)"), barcode)
  )
  tb_spill_prob <- metadata(sce)$spillover_est[[ch]]
  p_spill <- tb_bead %>%
    ggplot(aes(tfm(.data[[ch]]), color = .data$barcode)) +
    geom_density(adjust = 1, linewidth = 0.8) +
    geom_line(
      data = tb_spill_prob,
      aes(tfm(.data[[ch]]), .data$spill_prob),
      color = "black",
      linewidth = 0.8
    ) +
    ylab("density") +
    ggtitle("Beads Experiment")

  list(
    p_spill = p_spill,
    p_before_after = p_before_after
  )
}
