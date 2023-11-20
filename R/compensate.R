#' Compute spillover probability and correct for spillover
#'
#' @importFrom dplyr left_join filter mutate select bind_cols
#' @importFrom tibble tibble
#' @importFrom spatstat.geom ewcdf
#' @importFrom magrittr %>% %<>%
#' @importFrom tidyselect all_of
#' @importFrom stats binomial coef glm rbinom rpois runmed ecdf
#'
#' @param tb_real Data frame or tibble with proteins counts of real experiment
#' @param tb_bead Data frame or tibble with proteins counts of bead experiment
#' @param target_marker Marker name in real experiment
#' @param spillover_markers Marker names in bead experiment
#' @param runmed_k Integer width of median window for smoothing the ECDF
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
compensate <-
    function(tb_real,
             tb_bead,
             target_marker,
             spillover_markers,
             runmed_k, 
             n_iter = 1000) {
        # check if any beads
        tb_bead_keep <- tb_bead
        tb_bead <-
            tb_bead %>% filter(barcode != !!target_marker)
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
        tfm <- function(x)
            asinh(x / 5)
        smoothing <- function(pmf) {
            pmf_smooth <- runmed(pmf, k = runmed_k)
            pmf_smooth / sum(pmf_smooth)
        }
        epsilon <- 1/10^5
        all_markers <- c(target_marker, spillover_markers)
        y_min <- y_min_real
        y_max <- y_max_real
        
        # support for target marker
        tb_beads_pmf <- tibble(y = y_min:y_max)
        
        # collect pmf from beads
        for (marker in spillover_markers) {
            n <- nrow(filter(tb_bead, .data$barcode == marker))
            if (n > 0) {
                Fn <- tb_bead %>%
                    filter(.data$barcode == marker) %>%
                    pull(all_of(target_marker)) %>%
                    ecdf()
                tb <- tb_beads_pmf %>%
                    mutate(pmf = Fn(y) - Fn(y - 1)) %>%
                    mutate(pmf = smoothing(pmf)) %>%
                    select(y, pmf)
                names(tb) <- c("y", marker)
                tb_beads_pmf %<>% left_join(tb, by = "y")
            } else {
                tb <- tibble(y = y_min:y_max,
                             pmf = 1 / nrow(tb_beads_pmf))
                names(tb) <- c("y", marker)
                tb_beads_pmf %<>% left_join(tb, by = "y")
            }
        }
        
        # --------- step 1: initialize ---------
        
        # uniform prior probability
        n_markers <- length(all_markers)
        pi <- rep(NA, n_markers)
        pi[1] <- 0.9
        pi[2:n_markers] <- 0.1/(n_markers-1)
        names(pi) <- all_markers
        
        # add pmf from real cells
        Fn <- tb_real %>%
            pull(all_of(target_marker)) %>%
            ecdf()
        tb_real_pmf <- tibble(y = y_min:y_max) %>%
            mutate(pmf = Fn(y) - Fn(y - 1)) %>%
            mutate(pmf = smoothing(pmf)) %>%
            select(y, pmf)
        names(tb_real_pmf) <- c("y", target_marker)
        
        # join beads and real
        tb_pmf <- left_join(tb_beads_pmf, tb_real_pmf, by = "y")
        
        # --------- step 2: iterate ---------
        
        convergence <- matrix(nrow = n_iter, ncol = length(pi) + 1)
        colnames(convergence) <- c("iteration", names(pi))
        convergence[1, ] <- c(1, pi)
        
        for (i in seq(2, n_iter)) {
            # E-step
            
            # membership probabilities
            M <- tb_pmf %>% 
                select(all_of(all_markers)) %>% 
                as.matrix()
            PI <- matrix(rep(pi, each = nrow(M)), ncol = ncol(M))
            post_M <- PI * M
            post_M <- post_M / rowSums(post_M)
            # assign equal probability to counts without any signal
            post_M[is.nan(post_M)] <- 1 / ncol(post_M)
            
            # M-step
            
            # update prior probability
            y_pi <- bind_cols(y = tb_pmf$y, post_M)
            y_obsv <- tb_real %>% 
                select(y = all_of(target_marker))
            y_obsv <- left_join(y_obsv, y_pi, by = "y")
            pi <- y_obsv %>% select(-y)
            pi <- colSums(pi) / nrow(pi)
            
            # new weighted empirical density estimate
            Fn <- ewcdf(
                pull(y_obsv, y),
                weights = pull(y_obsv, all_of(target_marker))
                )
            tb_real_pmf <- tibble(y = y_min:y_max) %>%
                mutate(pmf = Fn(y) - Fn(y - 1)) %>%
                mutate(pmf = smoothing(pmf))
            names(tb_real_pmf) <- c("y", target_marker)
            
            # update join
            tb_pmf <- left_join(tb_beads_pmf, tb_real_pmf, by = "y")
            
            # keep track
            convergence[i, ] <- c(i, pi)
            
            prev <- convergence[i - 1, target_marker]
            curr <- convergence[i, target_marker]
            if (abs(prev - curr) < epsilon)
                break
        }
        
        convergence <- convergence[seq(i), ]
        
        # --------- spillover probability curve ---------
        
        # calculate posterior spillover probability for each cell
        M <- tb_pmf %>% 
            select(all_of(all_markers)) %>% 
            as.matrix()
        PI <- matrix(rep(pi, each = nrow(M)), ncol = ncol(M))
        post_M <- PI * M
        post_M <- post_M / rowSums(post_M)
        spill_prob <- 1 - post_M[, target_marker]
        tb_spill_prob <- tibble(y = y_min:y_max, spill_prob)
        names(tb_spill_prob) <- c(target_marker, "spill_prob")
        
        # compensate
        tb_compensate <- tb_real
        tb_compensate %<>% left_join(tb_spill_prob, by = target_marker)
        tb_compensate %<>% 
            mutate(spill_prob = ifelse(is.na(spill_prob), 0, spill_prob))
        tb_compensate %<>% mutate(spill = rbinom(
            n = nrow(tb_compensate),
            size = 1,
            prob = tb_compensate$spill_prob
        ))
        tb_compensate %<>%
            mutate(
                corrected = ifelse(
                    spill == 1, NA, .data[[target_marker]]))
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
