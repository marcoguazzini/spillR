load_beads_data <- function(){
  # loading the beads data
  # constants
  bc_key <- c(139, 141:156, 158:176)
left <- 2.5
right <- 2.7

# --------- experiment with beads ---------

sce_spill <- prepData(ss_exp)
sce_spill <- assignPrelim(sce_spill, bc_key, verbose = FALSE)
sce_spill <- applyCutoffs(estCutoffs(sce_spill))
sce_spill <- computeSpillmat(sce_spill)
counts_bead <- t(assay(sce_spill, "counts"))
counts_bead <- floor(counts_bead)
counts_bead <- as_tibble(counts_bead)

channel_names <- rowData(sce_spill)[,"channel_name"]
names(counts_bead) <- channel_names
counts_bead <- mutate(counts_bead, barcode = sce_spill$bc_id)
# --------- spillover matrix ---------
sm <- metadata(sce_spill)$spillover_matrix

list(counts_bead,sm)
}