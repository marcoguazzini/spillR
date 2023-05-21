# spillR: vignette

# Algorithm for real data.


# --------- experiment with beads ---------
bc_key <- c(139, 141:156, 158:176)
sce_bead <- prepData(ss_exp)
sce_bead <- assignPrelim(sce_bead, bc_key, verbose = FALSE)
sce_bead <- applyCutoffs(estCutoffs(sce_bead))
sce_bead <- computeSpillmat(sce_bead)

# --------- experiment with real cells ---------
data(mp_cells, package = "CATALYST")
sce <- prepData(mp_cells)

# --------- table for mapping markers and barcode ---------
marker_to_barc <- 
  rowData(sce_bead)[,c("channel_name", "is_bc")] %>%
  as_tibble %>%
  dplyr::filter(is_bc == TRUE) %>%
  mutate(barcode = bc_key) %>%
  select(marker = channel_name, barcode)

# --------- call compensate from compCytof package ---------
sce_spillr <- compCytof(sce, sce_bead, marker_to_barc, overwrite = FALSE)

# --------- 2d histogram from spillR ---------
 as <- c("counts", "exprs", "compcounts", "compexprs")
 chs <- c( "Yb171Di", "Yb173Di")
 ps <- lapply(as, function(a) 
     plotScatter(sce_spillr, chs, assay = a))
 plot_grid(plotlist = ps, nrow = 2)
 
 
```r 
# DAG representing the spillover
# set this for better image resolution fig.width=10, fig.height=7, out.width="100%"
spillR::spill_dag(spillR::load_spillover())
```





