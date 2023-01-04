plot_emit <- function(target) {
  barcode_target <- substr(target, 3, 5)
marker_selection <- names(smc[, target])[smc[, target] > 0]
barcodes <- sapply(marker_selection, function(marker) substr(marker, 3, 5))
bar <- barcodes[barcodes != barcode_target]
  emit_spill <-
    sapply(marker_selection, function(marker)
      substr(marker, 3, 5))
  
  y <- counts_spill %>%
    select_at(c(barc, "barcode")) %>%
    filter(barcode %in% emit_spill)
  colnames(y)[1] <- "x"
  # Plotting the various spillover distribution
  ggplot(y %>% filter(barcode != barc), aes(x)) +
    geom_histogram(bins = 50) +
    facet_wrap( ~ barcode)
}
