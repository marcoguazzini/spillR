load_real_data <- function(){
  # --------- experiment with real cells ---------
  data(mp_cells, package = "CATALYST")
  sce <- prepData(mp_cells)
  counts_real <- t(assay(sce, "counts"))
  counts_real <- floor(counts_real)
  colnames(counts_real) <- channel_names
  as_tibble(counts_real)
  
}