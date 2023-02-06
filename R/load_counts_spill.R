load_counts_spill <- function(){
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
  return(counts_spill)
}
