#' Compute spillover probability and correct for spillover
#'
#' @import dplyr
#' @import CATALYST
#' @importFrom magrittr %<>% %>%
#' @importFrom stats binomial coef glm rbinom rpois
#' @importFrom SummarizedExperiment assay rowData assay<-
#' @importFrom S4Vectors metadata
#' @export
#'
#' @param sce Single Cell Experiment for the real cells
#' @param sce_bead Single Cell Experiment for the bead experiment
#' @param marker_to_barc Table for barcodes in the beads experiment
#' @param overwrite logical; if TRUE data are overwritten if FALSE 
#'   data are saved in new columns
#'
#' @return A list of class \code{spillr} containing
#'   \item{tb_compensate}{corrected real cells}
#'   \item{tb_spill_prob}{probability curve}
#'   \item{target_marker}{input marker in real experiment}
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
#' as_tibble %>%
#' dplyr::filter(is_bc == TRUE) %>%
#' mutate(barcode = bc_key) %>%
#' dplyr::select(marker = channel_name, barcode)
#' spillR::compCytof(sce, sce_bead, marker_to_barc, overwrite = FALSE)
compCytof <- function(sce, sce_bead, marker_to_barc, overwrite = FALSE){
    if(!("marker" %in% colnames(marker_to_barc)))
        stop("marker_to_barc needs to have column marker")
    if(!("barcode" %in% colnames(marker_to_barc)))
        stop("marker_to_barc needs to have column barcode")
    # constants
    tfm <- function(x) asinh(x/5)
    # --------- experiment with beads ---------
    counts_bead <- t(SummarizedExperiment::assay(sce_bead, "counts"))
    counts_bead <- floor(counts_bead)
    counts_bead <- tibble::as_tibble(counts_bead)
    channel_names <- SummarizedExperiment::rowData(sce_bead)[,"channel_name"]
    names(counts_bead) <- channel_names
    counts_bead <-dplyr::mutate(counts_bead, barcode = sce_bead$bc_id)
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
                               dplyr::select(dplyr::all_of(target_marker)) %>%
                               dplyr::mutate(barcode = "none") %>%
                               dplyr::mutate(type = "real cells")
                           # rename barcodes to marker names
                           for(i in seq_len(length(spillover_barcodes))) {
                               ids <- tb_bead$barcode == spillover_barcodes[i]
                               tb_bead[ids, "barcode"] <- spillover_markers[i]
                               }
                           spillover_markers <- 
                               setdiff(spillover_markers, target_marker)
                           compensate(
                               tb_real, 
                               tb_bead, 
                               target_marker, 
                               spillover_markers)
                           }
                       )
    names(fit_list) <- rownames(sm)
    # saving the compensated data in a dataframe
    # adding the column for the markers not present in the beads experiment
    # in the beads experiment they do not attached three markers
    # "Ba138Di", "Ce140Di", "Gd157Di"
    # in general we can check the key barcode attached on the beads experiment
    # used channel in bead experiments
    used_channel <- sce_bead@rowRanges@elementMetadata@listData[["is_bc"]]
    c <- which(used_channel == FALSE)
    channels_out<- channel_names[c]
    # with the key we can check whether the
    data <- matrix(NA, nrow = nrow(counts_real), ncol = length(channel_names))
    data <- data.frame(data)
    for(i in seq_len(length(channel_names))){
        if(channel_names[i] %in% channels_out){
            data[,i] <- counts_real[,channel_names[i]]
            }
        else {
            data[,i] <- 
                unlist(fit_list[[channel_names[i]]]$tb_compensate[,"corrected"])
        }
        }
    colnames(data)<- channel_names
    # save compensated counts
    c <- ifelse(overwrite, "counts", "compcounts")
    SummarizedExperiment::assay(sce, c, FALSE) <- t(data)

    # save compensated transformed counts
    c <- ifelse(overwrite, "exprs", "compexprs")
    SummarizedExperiment::assay(sce, c, FALSE) <- t(tfm(data))
    return(sce)

}
