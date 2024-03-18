#' Compute spillover probability and correct for spillover from beads only
#'
#' @importFrom dplyr left_join filter mutate select bind_cols
#' @importFrom tibble tibble
#' @importFrom spatstat.geom ewcdf
#' @importFrom tidyselect all_of
#' @importFrom stats binomial coef glm rbinom rpois runmed ecdf
#'
#' @param tb_real Data frame or tibble with proteins counts of real experiment
#' @param tb_bead Data frame or tibble with proteins counts of bead experiment
#' @param target_marker Marker name in real experiment
#' @param spillover_markers Marker names in bead experiment
#' @param impute_value Value for counts that are declared as spillover
#'
#' @return A list of class \code{spillr} containing
#'   \item{tb_compensate}{corrected real cells}
#'   \item{tb_spill_prob}{probability curve}
#'   \item{convergence}{covergence table of EM algorithm}
#'   \item{tb_real}{input real cells}
#'   \item{tb_bead}{input bead cells}
#'   \item{target_marker}{input marker in real experiment}
#'   \item{spillover_markers}{input markers in bead experiment}
compensate_naive <- function(tb_real, tb_bead, target_marker, spillover_markers,
                             impute_value = NA) {
    # check if any beads
    if (nrow(tb_bead) == 0) {
        warning("no beads")
        return(NULL)
    }

    # check if real and bead overlap
    y_min_real <- min(tb_real[, target_marker])
    y_max_real <- max(tb_real[, target_marker])
    y_min_bead <- min(tb_bead[, target_marker])
    y_max_bead <- max(tb_bead[, target_marker])
    set <- intersect(y_min_real:y_max_real, y_min_bead:y_max_bead)
    if (length(set) == 0) {
        warning("no overlap between real cells and beads")
        return(NULL)
    }

    # parameters and helper functions
    epsilon <- 1 / 10^15
    all_markers <- c(target_marker, spillover_markers)
    y_min <- y_min_real
    y_max <- y_max_real

    # support for target marker
    tb_beads_pmf <- tibble(y = y_min:y_max)

    # collect pmf from beads
    tb_beads_smooth <- lapply(all_markers, function(marker) {
        
        n <- nrow(filter(tb_bead, barcode == marker))
        tb <- tibble(
            y = y_min:y_max,
            pmf = 1 / nrow(tb_beads_pmf)
        )
        if (n > 0) {
            
            y <- tb_bead |>
                filter(barcode == marker) |>
                pull(all_of(target_marker))
            fit <- density(y, from = y_min, to = y_max)
            f <- approxfun(fit$x, fit$y)
            
            tb <- tb_beads_pmf |>
                mutate(pmf = f(y)) |>
                mutate(pmf = pmf/sum(pmf)) |>
                select(y, pmf)
        }
        names(tb)[2] <- marker
        tb[, marker]
        
    }) |> bind_cols()

    tb_beads_pmf <- bind_cols(tb_beads_pmf, tb_beads_smooth)

    # --------- spillover probability curve ---------

    # calculate posterior spillover probability for each cell
    M <- tb_beads_pmf |>
        select(all_of(all_markers)) |>
        as.matrix()
    M_norm <- sweep(M, 1, rowSums(M), FUN = "/")
    spill_prob <- 1 - M_norm[, target_marker]
    no_signal <- apply(M, 1, max) < epsilon
    spill_prob[no_signal] <- 0
    
    tb_spill_prob <- tibble(y = y_min:y_max, spill_prob)
    names(tb_spill_prob) <- c(target_marker, "spill_prob")

    # compensate
    tb_compensate <- tb_real
    tb_compensate <- left_join(tb_compensate, tb_spill_prob, by = target_marker)
    tb_compensate <- mutate(
        tb_compensate,
        spill_prob = ifelse(is.na(spill_prob), 0, spill_prob)
    )
    tb_compensate <- mutate(
        tb_compensate,
        spill = rbinom(
            n = nrow(tb_compensate),
            size = 1,
            prob = tb_compensate$spill_prob
        )
    )
    tb_compensate <- mutate(
        tb_compensate,
        corrected = ifelse(spill == 1, impute_value, .data[[target_marker]])
    )
    names(tb_compensate)[1] <- "uncorrected"

    # return spillr object
    res <- NULL
    res$tb_compensate <- tb_compensate
    res$tb_spill_prob <- tb_spill_prob
    res$tb_real <- tb_real
    res$tb_bead <- tb_bead
    res$target_marker <- target_marker
    res$spillover_markers <- spillover_markers
    class(res) <- "spillr"
    res
}
