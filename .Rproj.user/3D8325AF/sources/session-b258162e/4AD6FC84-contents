#' Compute spillover probability and correct for spillover
#'
#' @import magrittr
#' @import dplyr
#' @export
#'
#' @param tb_real Data frame or tibble with proteins counts of real experiment
#' @param tb_bead Data frame or tibble with proteins counts of bead experiment
#' @param target_marker Marker name in real experiment
#' @param spillover_markers Marker names in bead experiment
#'
#' @return A list of class \code{spillr} containing
#'   \item{tb_compensate}{corrected real cells}
#'   \item{tb_spill_prob}{probability curve}
#'   \item{convergence}{covergence table of EM algorithm}
#'   \item{tb_real}{input real cells}
#'   \item{tb_bead}{input bead cells}
#'   \item{target_marker}{input marker in real experiment}
#'   \item{spillover_markers}{input markers in bead experiment}
#'
#' @examples
#' set.seed(23)
#' tb_real <- generate_real()
#' tb_bead <- generate_bead()
#' target_marker <- "A"
#' spillover_marker <- "B"
#' spillR::compensate(tb_real, tb_bead, "A", "B")
compensate <- function(tb_real, tb_bead, target_marker, spillover_markers) {
  
  # --------- mixture method ---------
  
  # option 1: fit polynomial poisson regression
  # degree_bead <- 4
  # degree_real <- 8
  
  # option 2: smoother
  # degree_bead <- 11
  # degree_real <- 51
  
  # option 3: local regression
  # degree_bead <- 0.01
  # degree_real <- 0.05
  
  # option 4: no smoothing
  # no hyperparameters
  
  tfm <- function(x) asinh(x/5)
  all_markers <- c(target_marker, spillover_markers)
  y_min <- tb_real %>% pull(all_of(target_marker)) %>% min()
  y_max <- tb_real %>% pull(all_of(target_marker)) %>% max()
  
  denoise <- function(y, y_min = min(y), y_max = max(y)) {
    
    # frequency table
    tb_obsv <- tibble(y) %>% group_by(y) %>% tally()
    y_min_obsv <- min(tb_obsv$y)
    y_max_obsv <- max(tb_obsv$y)
    tb_pred <- tibble(y = y_min:y_max)
    tb_pred %<>% left_join(tb_obsv, by = "y")
    
    # padding with zero outside of data support
    tb_pred %<>% mutate(n = ifelse(y > y_max_obsv | y < y_min_obsv, 0, n))
    
    # option 1
    # fit <- glm(n ~ poly(y, degree = degree, raw = TRUE),
    #            data = tb_pred,
    #            family = "poisson")
    # tb_pred %<>% mutate(lambda = exp(predict(fit, newdata = tb_pred)))
    
    # option 2
    # n_smooth <- runmed(tb_pred$n, k = degree)
    # n_smooth[is.na(n_smooth)] <- 0
    # tb_pred %<>% mutate(lambda = n_smooth)
    
    # option 3
    # fit <- loess(n ~ y, data = tb_pred, span = degree)
    # tb_pred %<>% mutate(lambda = predict(fit, newdata = tb_pred))
    # tb_pred %<>% mutate(lambda = ifelse(lambda < 0, 0, lambda))
    
    # option 4
    tb_pred$n[is.na(tb_pred$n)] <- 0
    tb_pred %<>% mutate(lambda = n)
    
    # normalize
    tb_pred %>% mutate(pmf = lambda/sum(tb_pred$lambda))
    
  }
  
  # support for target marker
  tb_beads_pmf <- tibble(y = y_min:y_max)
  
  # collect pmf from beads
  for(marker in spillover_markers) {
    
    tb <- tb_bead %>% 
      dplyr::filter(barcode == marker) %>% 
      pull(all_of(target_marker)) %>% 
      denoise(y_min = y_min, y_max = y_max) %>% 
      dplyr::select(y, pmf)
    names(tb) <- c("y", marker)
    tb_beads_pmf %<>% left_join(tb, by = "y")
    
  }
  
  # --------- step 1: initialize ---------
  
  # prior probability
  pi <- rep(1, length(all_markers))
  pi <- pi/length(pi)
  names(pi) <- all_markers
  
  # add pmf from real cells
  tb_real_pmf <- tb_real %>% 
    pull(all_of(target_marker)) %>% 
    denoise(y_min = y_min, y_max = y_max) %>% 
    dplyr::select(y, pmf)
  names(tb_real_pmf) <- c("y", target_marker)
  
  # join beads and real
  tb_pmf <- left_join(tb_beads_pmf, tb_real_pmf, by = "y")
  
  # --------- step 2: iterate ---------
  
  n_iter <- 10
  convergence <- matrix(nrow = n_iter, ncol = length(pi)+1)
  colnames(convergence) <- c("iteration", names(pi))
  
  for(i in 1:n_iter) {
    # remove counts without any signal
    total <- rowSums(dplyr::select(tb_pmf, all_of(all_markers)))
    tb_pmf %<>% mutate(total = total)
    tb_pmf %<>% dplyr::filter(total > 0)
    tb_pmf %<>% dplyr::select(-total)
    
    # E-step
    
    # membership probabilities
    M <- tb_pmf %>% dplyr::select(all_of(all_markers)) %>% as.matrix()
    PI <- matrix(rep(pi, each = nrow(M)), ncol = ncol(M))
    post_M <- PI * M
    post_M <- post_M/rowSums(post_M)
    
    # M-step
    
    # update prior probability
    pi <- colSums(post_M)/nrow(post_M)
    
    # stochastic EM: 
    # assigns each observation to a class with the highest posterior probability
    #class <- apply(post_M, 1, function(Mi) sample(colnames(post_M), size = 1, prob = Mi))
    # categorical EM:
    # assigns each observation randomly based on posterior probabilities
    class <- all_markers[apply(post_M, 1, which.max)]
    
    # update pmf from real cells
    if(sum(class == target_marker) > 0) {
      
      ys <- tb_pmf[class == target_marker, ] %>% pull(y)
      tb_real_pmf <- tb_real %>% 
        dplyr::filter(.data[[target_marker]] %in% ys) %>%
        pull(all_of(target_marker)) %>% 
        denoise(y_min = y_min, y_max = y_max) %>% 
        dplyr::select(y, pmf)
      
    } else {
      
      # if no signal, then use uniform distribution  
      tb_real_pmf <- tb_beads_pmf %>% 
        dplyr::select(y) %>% 
        mutate(pmf = 1/nrow(tb_beads_pmf))
      
    }
    names(tb_real_pmf) <- c("y", target_marker)
    
    # update join
    tb_pmf <- left_join(tb_beads_pmf, tb_real_pmf, by = "y")
    
    # keep track
    convergence[i,] <- c(i, pi)
  }
  
  # --------- spillover probability curve ---------
  
  # calculate posterior spillover probability for each cell
  M <- tb_pmf %>% dplyr::select(all_of(all_markers)) %>% as.matrix()
  PI <- matrix(rep(pi, each = nrow(M)), ncol = ncol(M))
  post_M <- PI * M
  post_M <- post_M/rowSums(post_M)
  spill_prob <- 1-post_M[,target_marker]
  tb_spill_prob <- dplyr::select(tb_pmf, y) %>% mutate(spill_prob = spill_prob)
  tb_spill_prob %<>% mutate(spill_prob = if_else(is.na(spill_prob), 0, spill_prob))
  names(tb_spill_prob) <- c(target_marker, "spill_prob")
  
  # compensate
  tb_compensate <- tb_real
  tb_compensate %<>% left_join(tb_spill_prob, by = target_marker)
  tb_compensate %<>% mutate(
    spill = rbinom(n = nrow(tb_compensate), 
                   size = 1, 
                   prob = tb_compensate$spill_prob)
  )
  tb_compensate %<>% mutate(corrected = ifelse(spill == 1, NA, .data[[target_marker]]))
  
  names(tb_compensate)[1] = "uncorrected"
  
  # postprocess spillover probabilities
  tb_spill_prob %<>% mutate(y_tfm = tfm(.data[[target_marker]]))
  fit <- glm(spill_prob ~ y_tfm, family = binomial, data = tb_spill_prob)
  inverse_logit <- function(fit, x) {
    hat <- coef(fit)[1] + coef(fit)[2]*x
    1/(1+exp(-hat))
  }
  tb_spill_prob %<>% mutate(spill_prob_smooth = inverse_logit(fit, y_tfm))
  
  # return spillr object
  res <- NULL
  res$tb_compensate <- tb_compensate
  res$tb_spill_prob <- tb_spill_prob
  res$convergence <- convergence
  res$tb_real <- tb_real
  res$tb_bead <- tb_bead
  res$target_marker <- target_marker
  res$spillover_markers <- spillover_markers
  class(res) <- "spillr"
  res
  
}
