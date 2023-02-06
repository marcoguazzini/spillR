plot_emit <- function(target) {
  barcode_target <- substr(target, 3, 5)
marker_selection <- names(smc[, target])[smc[, target] > 0]
barcodes <- sapply(marker_selection, function(marker) substr(marker, 3, 5))
bar <- barcodes[barcodes != barcode_target]
  emit_spill <-
    sapply(marker_selection, function(marker)
      substr(marker, 3, 5))
  data(ss_exp, package = "CATALYST")
# specify mass channels stained for & debarcode
bc_ms <- c(139, 141:156, 158:176)
sce_spill <- prepData(ss_exp)
sce_spill <- assignPrelim(sce_spill, bc_ms, verbose = FALSE)
sce_spill <- applyCutoffs(estCutoffs(sce_spill))
counts_spill <- t(assay(sce_spill, "counts"))
counts_spill <- floor(counts_spill)
channel_names_spill <- as.matrix(rowData(sce_spill))[,"channel_name"]
colnames(counts_spill) <- channel_names_spill
# There are a total of 36 different experiments.
# Each experiment has a barcode to identify it. 
# In each experiment only one marker is attached, 
# the others are due to spillover.
counts_spill %<>% as_tibble()
names(counts_spill) <- sapply(names(counts_spill), 
                              function(marker) substr(marker, 3, 5))
counts_spill %<>% mutate(barcode = sce_spill$bc_id)
  y <- counts_spill %>%
    select_at(c(barcode_target, "barcode")) %>%
    filter(barcode %in% emit_spill)
  colnames(y)[1] <- "x"
  # Plotting the various spillover distribution
  ggplot(y %>% filter(barcode != barcode_target), aes(x)) +
    geom_histogram(bins = 50) +
    facet_wrap( ~ barcode)
}
