#' Compute spillover probability and correct for spillover
#'
#' @importFrom CATALYST prepData assignPrelim applyCutoffs computeSpillmat
#'                      prepData
#' @importFrom dplyr as_tibble pull mutate select
#' @importFrom stats binomial coef glm rbinom rpois
#' @importFrom SummarizedExperiment assay rowData assay<-
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom parallel mclapply
#' @export
#'
#' @param sce \code{\link[SingleCellExperiment]{SingleCellExperiment}} for the
#'   real cells
#' @param sce_bead \code{\link[SingleCellExperiment]{SingleCellExperiment}} for
#'   the bead experiment
#' @param marker_to_barc Table that maps the marker to the barcode
#'   in the beads experiment
#' @param impute_value Imputed value for counts that are declared as spillover
#' @param overwrite logical; if TRUE data are overwritten if FALSE
#'     data are saved in new columns
#' @param n_cores Number of computing cores
#' @param naive logical; if TRUE use the naive version
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'
#' @examples
#' library(CATALYST)
#' library(dplyr)
#' bc_key <- c(139, 141:156, 158:176)
#' sce_bead <- prepData(ss_exp)
#' sce_bead <- assignPrelim(sce_bead, bc_key, verbose = FALSE)
#' sce_bead <- applyCutoffs(estCutoffs(sce_bead))
#' sce_bead <- computeSpillmat(sce_bead)
#' data(mp_cells, package = "CATALYST")
#' sce <- prepData(mp_cells)
#' marker_to_barc <- rowData(sce_bead)[, c("channel_name", "is_bc")] |>
#'     as_tibble() |>
#'     filter(is_bc == TRUE) |>
#'     mutate(barcode = bc_key) |>
#'     select(marker = channel_name, barcode)
#' spillR::compCytof(sce, sce_bead, marker_to_barc, impute_value = NA)
compCytof <- function(sce, sce_bead, marker_to_barc, impute_value, 
                      overwrite = FALSE, n_cores = 1, naive = FALSE) {
    if (!("marker" %in% colnames(marker_to_barc))) {
        stop("marker_to_barc needs to have column marker")
    }
    if (!("barcode" %in% colnames(marker_to_barc))) {
        stop("marker_to_barc needs to have column barcode")
    }

    # --------- experiment with beads ---------

    counts_bead <- t(assay(sce_bead, "counts"))
    counts_bead <- floor(counts_bead)
    counts_bead <- as_tibble(counts_bead)

    channel_names <- rowData(sce_bead)[, "channel_name"]
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

    fit_list <- mclapply(
        rownames(sm),
        function(target_marker) {
            spillover_markers <- names(which(sm[, target_marker] > 0))
            spillover_barcodes <- marker_to_barc |>
                filter(marker %in% spillover_markers) |>
                select("barcode") |>
                pull()

            tb_bead <- counts_bead |>
                filter(barcode %in% spillover_barcodes) |>
                select(all_of(c(target_marker, "barcode"))) |>
                mutate(type = "beads")

            tb_real <- counts_real |>
                select(all_of(target_marker)) |>
                mutate(barcode = "none") |>
                mutate(type = "real cells")

            # rename barcodes to marker names
            ids <- vapply(
                tb_bead$barcode,
                function(bc) which(bc == spillover_barcodes), numeric(1)
            )
            tb_bead$barcode <- spillover_markers[ids]

            spillover_markers <- setdiff(spillover_markers, target_marker)
            
            if(!naive) {
                compensate(tb_real, tb_bead, target_marker, 
                           spillover_markers, impute_value)
            } else {
                compensate_naive(tb_real, tb_bead, target_marker, 
                                 spillover_markers, impute_value)
            }
        }, 
        mc.cores = n_cores
    )
    names(fit_list) <- rownames(sm)

    # --------- save results in SingleCellExperiment class ---------

    # find inactive channels
    used_channel <- rowData(sce_bead)$is_bc
    channels_out <- channel_names[used_channel == FALSE]
    channels_null <- names(which(vapply(fit_list, is.null, logical(1))))
    channels_out <- union(channels_out, channels_null)
    channels_in <- setdiff(channel_names, channels_out)

    # prepare new assay matrices
    data <- matrix(NA, nrow = nrow(counts_real), ncol = length(channel_names))
    data <- data.frame(data)
    colnames(data) <- channel_names
    spillprob <- data

    # no correction
    data[, channels_out] <- counts_real[, channels_out]

    # keep the corrected counts
    data[, channels_in] <- lapply(
        fit_list[channels_in],
        function(fit) fit$tb_compensate$corrected
    )

    # keep the smoothed spillover probability for diagnostic plots
    spillprob[, channels_in] <- lapply(
        fit_list[channels_in],
        function(fit) fit$tb_compensate$spill_prob
    )

    # save compensated counts
    c <- ifelse(overwrite, "counts", "compcounts")
    assay(sce, c, FALSE) <- t(data)

    # save compensated transformed counts
    c <- ifelse(overwrite, "exprs", "compexprs")
    assay(sce, c, FALSE) <- t(tfm(data))

    # save spillover probabilities
    assay(sce, "spillprob", FALSE) <- t(spillprob)

    # add spillover meta data for diagnostic plots
    beads_distr <- lapply(fit_list, function(fit) fit$tb_bead)
    spillover_est <- lapply(fit_list, function(fit) fit$tb_spill_prob)
    metadata(sce)$beads_distr <- beads_distr
    metadata(sce)$spillover_est <- spillover_est

    return(sce)
}
