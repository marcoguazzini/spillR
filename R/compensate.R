#' Compute spillover probability and correct for spillover
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
#' @param n_iter Maximum number of EM steps
#'
#' @return A list of class \code{spillr} containing
#'   \item{tb_compensate}{corrected real cells}
#'   \item{tb_spill_prob}{probability curve}
#'   \item{convergence}{covergence table of EM algorithm}
#'   \item{tb_real}{input real cells}
#'   \item{tb_bead}{input bead cells}
#'   \item{target_marker}{input marker in real experiment}
#'   \item{spillover_markers}{input markers in bead experiment}
compensate <- function(tb_real, tb_bead, target_marker, spillover_markers,
                       impute_value = NA, n_iter = 1000) {
    # check if any beads
    tb_bead_keep <- tb_bead
    tb_bead <-
        tb_bead |> filter(barcode != !!target_marker)
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
    epsilon <- 1 / 10^5
    all_markers <- c(target_marker, spillover_markers)
    y_min <- y_min_real
    y_max <- y_max_real

    # support for target marker
    tb_beads_pmf <- tibble(y = y_min:y_max)

    # collect pmf from beads
    tb_beads_smooth <- lapply(spillover_markers, function(marker) {
        
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

    # --------- step 1: initialize ---------

    # uniform prior probability
    n_markers <- length(all_markers)
    pi <- rep(NA, n_markers)
    pi[1] <- 0.9
    pi[2:n_markers] <- 0.1 / (n_markers - 1)
    names(pi) <- all_markers

    # add pmf from real cells
    y <- pull(tb_real, all_of(target_marker))
    fit <- density(y, from = y_min, to = y_max)
    f <- approxfun(fit$x, fit$y)
    
    tb_real_pmf <- tibble(y = y_min:y_max) |>
        mutate(pmf = f(y)) |>
        mutate(pmf = pmf/sum(pmf)) |>
        select(y, pmf)
    names(tb_real_pmf) <- c("y", target_marker)

    # join beads and real
    tb_pmf <- left_join(tb_beads_pmf, tb_real_pmf, by = "y")

    # --------- step 2: iterate ---------

    convergence <- matrix(nrow = n_iter, ncol = length(pi) + 1)
    colnames(convergence) <- c("iteration", names(pi))
    convergence[1, ] <- c(1, pi)

    prev <- epsilon
    curr <- -epsilon
    i <- 2
    while (i <= n_iter & abs(prev - curr) > epsilon) {
        # E-step

        # membership probabilities
        M <- tb_pmf |>
            select(all_of(all_markers)) |>
            as.matrix()
        PI <- matrix(rep(pi, each = nrow(M)), ncol = ncol(M))
        post_M <- PI * M
        post_M <- post_M / rowSums(post_M)
        # assign equal probability to counts without any signal
        post_M[is.nan(post_M)] <- 1 / ncol(post_M)

        # M-step

        # update prior probability
        y_pi <- bind_cols(y = tb_pmf$y, post_M)
        y_obsv <- tb_real |>
            select(y = all_of(target_marker))
        y_obsv <- left_join(y_obsv, y_pi, by = "y")
        pi <- y_obsv |> select(-y)
        pi <- colSums(pi) / nrow(pi)

        # new weighted empirical density estimate
        y <- pull(y_obsv, y)
        weights = pull(y_obsv, all_of(target_marker))
        fit <- density(y, from = y_min, to = y_max, weights = weights)
        f <- approxfun(fit$x, fit$y)
        
        tb_real_pmf <- tibble(y = y_min:y_max) |>
            mutate(pmf = f(y)) |>
            mutate(pmf = pmf/sum(pmf)) |>
            select(y, pmf)
        names(tb_real_pmf) <- c("y", target_marker)

        # update join
        tb_pmf <- left_join(tb_beads_pmf, tb_real_pmf, by = "y")

        # keep track
        convergence[i, ] <- c(i, pi)
        prev <- convergence[i - 1, target_marker]
        curr <- convergence[i, target_marker]
        i <- i + 1
    }

    convergence <- convergence[seq(i - 1), ]

    # --------- spillover probability curve ---------

    # calculate posterior spillover probability for each cell
    M <- tb_pmf |>
        select(all_of(all_markers)) |>
        as.matrix()
    PI <- matrix(rep(pi, each = nrow(M)), ncol = ncol(M))
    post_M <- PI * M
    post_M <- post_M / rowSums(post_M)
    spill_prob <- 1 - post_M[, target_marker]
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
    res$convergence <- convergence
    res$tb_real <- tb_real
    res$tb_bead <- tb_bead_keep
    res$target_marker <- target_marker
    res$spillover_markers <- spillover_markers
    class(res) <- "spillr"
    res
}
