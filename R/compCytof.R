#' Compute spillover probability and correct for spillover
#'
#' @import dplyr
#' @import CATALYST
#' @importFrom magrittr %<>% %>%
#' @importFrom stats binomial coef glm rbinom rpois
#' @importFrom SummarizedExperiment assay rowData assay<-
#' @importFrom S4Vectors metadata metadata<-
#' @export
#'
#' @param sce Single Cell Experiment for the real cells
#' @param sce_bead Single Cell Experiment for the bead experiment
#' @param marker_to_barc Table for barcodes in the beads experiment
#' @param overwrite logical; if TRUE data are overwritten if FALSE 
#'     data are saved in new columns
#'
#' @return Compensates the input \code{\link{SummarizedExperiment}} or, 
#'     if \code{x} is a character string, 
#'     all FCS files in the specified location.
#'
#' @examples
#' library(CATALYST)
#' library(dplyr)
#' library(magrittr)
#' bc_key <- c(139, 141:156, 158:176)
#' sce_bead <- prepData(ss_exp)
#' sce_bead <- assignPrelim(sce_bead, bc_key, verbose = FALSE)
#' sce_bead <- applyCutoffs(estCutoffs(sce_bead))
#' sce_bead <- computeSpillmat(sce_bead)
#' data(mp_cells, package = "CATALYST")
#' sce <- prepData(mp_cells)
#' marker_to_barc <- rowData(sce_bead)[,c("channel_name", "is_bc")] %>%
#'     as_tibble %>%
#'     dplyr::filter(is_bc == TRUE) %>%
#'     mutate(barcode = bc_key) %>%
#'     dplyr::select(marker = channel_name, barcode)
#' spillR::compCytof(sce, sce_bead, marker_to_barc, overwrite = FALSE)
compCytof <- function(sce, sce_bead, marker_to_barc, overwrite = FALSE){
    if(!("marker" %in% colnames(marker_to_barc)))
        stop("marker_to_barc needs to have column marker")
    if(!("barcode" %in% colnames(marker_to_barc)))
        stop("marker_to_barc needs to have column barcode")
    tfm <- function(x) asinh(x/5) 
    # --------- experiment with beads ---------
    counts_bead <- t(SummarizedExperiment::assay(sce_bead, "counts"))
    counts_bead <- floor(counts_bead)
    counts_bead <- tibble::as_tibble(counts_bead)
    channel_names <- SummarizedExperiment::rowData(sce_bead)[,"channel_name"]
    names(counts_bead) <- channel_names
    counts_bead <- dplyr::mutate(counts_bead, barcode = sce_bead$bc_id)
    # what markers spillover into target marker?
    sm <- S4Vectors::metadata(sce_bead)$spillover_matrix
    
    # --------- experiment with real cells ---------
    
    channel_names <- SummarizedExperiment::rowData(sce)[, "channel_name"]
    counts_real <- t(SummarizedExperiment::assay(sce, "counts"))
    counts_real <- floor(counts_real)
    colnames(counts_real) <- channel_names
    counts_real <- tibble::as_tibble(counts_real)
    
    # --------- iterator over markers ---------
    fit_list <- lapply(rownames(sm), 
                       function(target_marker) {
                           spillover_markers <- 
                               names(which(sm[,target_marker] > 0))
                           spillover_barcodes <- marker_to_barc %>% 
                               dplyr::filter(
                                   .data$marker %in% spillover_markers) %>%
                               dplyr::select("barcode")%>%
                               dplyr::pull()
                           tb_bead <- counts_bead %>% 
                               dplyr::filter(
                                   .data$barcode %in% spillover_barcodes) %>% 
                               dplyr::select(
                                   dplyr::all_of(
                                       c(target_marker, "barcode"))) %>% 
                               dplyr::mutate(type = "beads")
                           tb_real <- counts_real %>% 
                               dplyr::select(all_of(target_marker)) %>%
                               dplyr::mutate(barcode = "none") %>%
                               dplyr::mutate(type = "real cells")
                           
                           # rename barcodes to marker names
                           for(i in seq_len(length(spillover_barcodes))) {
                               ids <- tb_bead$barcode == spillover_barcodes[i]
                               tb_bead[ids, "barcode"] <- spillover_markers[i]
                           }
                           spillover_markers <- 
                               setdiff(spillover_markers, target_marker)
                           compensate(tb_real, 
                                      tb_bead, 
                                      target_marker, 
                                      spillover_markers)
                       }
    )
    names(fit_list) <- rownames(sm)
    
    # --------- save results in SingleCellExperiment class ---------
    
    # find inactive channels
    used_channel <- SummarizedExperiment::rowData(sce_bead)$is_bc
    channels_out <- channel_names[used_channel == FALSE]
    
    # prepare new assay matrices
    data <- matrix(NA, nrow=nrow(counts_real), ncol=length(channel_names))
    data <- data.frame(data)
    colnames(data) <- channel_names
    spillprob <- data
    rep(NA, ncol(data))
    for(i in seq_len(length(channel_names))){
        if(channel_names[i] %in% channels_out){
            # no correct
            data[,i] <- counts_real[,channel_names[i]]
        }
        else
        {
            # keep the corrected counts
            tb_compensate <- fit_list[[channel_names[i]]]$tb_compensate
            data[,i]      <- tb_compensate$corrected
            
            # keep the smoothed spillover probability for diagnostic plots
            spillprob[,i] <- tb_compensate$spill_prob
        }
    }
    
    # save compensated counts
    c <- ifelse(overwrite, "counts", "compcounts")
    SummarizedExperiment::assay(sce, c, FALSE) <- t(data)
    
    # save compensated transformed counts
    c <- ifelse(overwrite, "exprs", "compexprs")
    SummarizedExperiment::assay(sce, c, FALSE) <- t(tfm(data))
    
    # save spillover probabilities
    SummarizedExperiment::assay(sce, "spillprob", FALSE) <- t(spillprob)
    
    # add spillover meta data for diagnostic plots
    beads_distr <- lapply(fit_list, function(fit) fit$tb_bead)
    spillover_est <- lapply(fit_list, function(fit) fit$tb_spill_prob)
    S4Vectors::metadata(sce)$beads_distr <- beads_distr
    S4Vectors::metadata(sce)$spillover_est <- spillover_est
    
    return(sce)
    
}
