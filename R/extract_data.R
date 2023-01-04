extract_data <- function()
{
  data(mp_cells, package = "CATALYST")
  sce <- prepData(mp_cells)
  channel_names <- as.matrix(rowData(sce))[,"channel_name"]
  counts <- t(assay(sce, "counts"))
  counts <- floor(counts)
  colnames(counts) <- channel_names
  counts
}