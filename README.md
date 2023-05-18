# spillR

 # Algorithm for real data.
 
How to use our package on real data:
  
 ```{r spillr-vignette, echo = FALSE, warning = FALSE, message = FALSE}
# constants
bc_key <- c(139, 141:156, 158:176)

# --------- experiment with beads ---------

sce_bead <- prepData(ss_exp)
sce_bead <- assignPrelim(sce_bead, bc_key, verbose = FALSE)
sce_bead <- applyCutoffs(estCutoffs(sce_bead))
sce_bead <- computeSpillmat(sce_bead)

# --------- experiment with real cells ---------

data(mp_cells, package = "CATALYST")
sce <- prepData(mp_cells)

# --------- call compensate from spillR compCytof package ---------
sce <- spillR::compCytof(sce, sce_bead, overwrite = FALSE) 
```

 
 
 
 
```r 
# DAG representing the spillover
# set this for better image resolution fig.width=10, fig.height=7, out.width="100%"
spillR::spill_dag(spillR::load_spillover())
```





