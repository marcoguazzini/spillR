#' Compute spillover probability and correct for spillover
#'
#' @import magrittr
#' @import dplyr
#' @import CATALYST
#' @import flowCore
#' @import tidyr
#' @export
#'
#' @param sce Single Cell Experiment for the real cells
#' @param sce_bead Single Cell Experiment for the bead experiment
#'
#' @return A list of class \code{spillr} containing
#'   \item{tb_compensate}{corrected real cells}
#'   \item{tb_spill_prob}{probability curve}
#'   \item{target_marker}{input marker in real experiment}
#'
#' @examples
#' sce <- prepData(mp_cells)    # real cells
#' sce_bead <- prepData(ss_exp) # beads
#' sce <- spillR::compCytof(sce, sce_bead, overwrite = FALSE)
compCytof <- function(sce, sce_bead, overwrite = FALSE){
  
  # modified from methods.Rmd
  
  # --------- experiment with beads ---------
  
  counts_bead <- t(assay(sce_bead, "counts"))
  counts_bead <- floor(counts_bead)
  counts_bead <- as_tibble(counts_bead)
  
  channel_names <- rowData(sce_bead)[,"channel_name"]
  names(counts_bead) <- channel_names
  counts_bead <- mutate(counts_bead, barcode = sce_bead$bc_id)
  
  # what markers spillover into target marker?
  sm <- metadata(sce_bead)$spillover_matrix
  
  # --------- experiment with real cells ---------
  
  channel_names <- rowData(sce)[, "channel_name"]
  counts_real <- t(assay(sce, "counts"))
  counts_real <- floor(counts_real)
  colnames(counts_real) <- channel_names
  counts_real <- as_tibble(counts_real)
  
  # --------- iterator over markers ---------
  
  # constants
  tfm <- function(x) asinh(x/5)
  marker_to_code <- function(marker) as.integer(substr(marker, 3, 5)) # TODO: this won't work in general
  
  fit_list <- lapply(rownames(sm), 
    function(target_marker) {
      
      spillover_markers <- names(which(sm[,target_marker] > 0))
      spillover_barcodes <- marker_to_code(spillover_markers)
      tb_bead <- counts_bead %>% 
        filter(barcode %in% spillover_barcodes) %>% 
        select(all_of(c(target_marker, "barcode"))) %>% 
        mutate(type = "beads")
      
      tb_real <- counts_real %>% 
        select(all_of(target_marker)) %>%
        mutate(barcode = "none") %>%
        mutate(type = "real cells")
      
      # rename barcodes to marker names
      for(i in seq(length(spillover_barcodes))) {
        ids <- tb_bead$barcode == spillover_barcodes[i]
        tb_bead[ids, "barcode"] <- spillover_markers[i]
      }
      
      spillover_markers <- setdiff(spillover_markers, target_marker)
      compensate(tb_real, tb_bead, target_marker, spillover_markers)
      
      }
    )
  
  names(fit_list) <- rownames(sm)
  fit_list 
  # saving the compensated data in a dataframe
  # adding the column for the markers not present in the beads experiment
  # in the beads experiment they do not attached three markers "Ba138Di", "Ce140Di", "Gd157Di"
  # in general we can check the key barcode attached on the beads experiment (bc_key)
  bc_key <- names(sce_bead@metadata[["bc_key"]])
  bc_real <- unlist(lapply(channel_names, function(channel) marker_to_code(channel)))
  bc_diff <- setdiff(bc_real,bc_key)
  channels_out <- rep(NA,length(bc_diff))
  for(i in 1:length(bc_diff)){
    c <- which(bc_real == bc_diff[i])
    channels_out[i]<- channel_names[c]
  }
  # with the key we can check whether the 
  data <- matrix(NA, nrow = nrow(counts_real), ncol = length(channel_names))
  data <- data.frame(data)
  for(i in 1:length(channel_names)){
    if(channel_names[i] %in% channels_out){
      data[,i] <- counts_real[,channel_names[i]]
    }
    else
    {
      data[,i] <- unlist(fit_list[[channel_names[i]]]$tb_compensate[,"corrected"])
    }
  }
  colnames(data)<- channel_names
  
  # save compensated counts
  c <- ifelse(overwrite, assay, "compcounts")
  assay(sce, c, FALSE) <- t(data)
  
  # save compensated transformed counts
  c <- ifelse(overwrite, "exprs", "compexprs")
  assay(sce, c, FALSE) <- t(tfm(data))
  return(sce)

}
