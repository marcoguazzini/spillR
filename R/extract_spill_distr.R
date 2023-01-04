extract_spill_distr <- function(){
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
  counts_spill
}