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
  
  # helper functions
  # helper functions
  tfm <- function(x) asinh(x/5)
  colorscale = scale_fill_gradientn(
    colors = rev(brewer.pal(9, "YlGnBu")),
    values = c(0, exp(seq(-5, 0, length.out = 100)))
  )
  marker_to_code <- function(marker) as.integer(substr(marker, 3, 5))
  
  # using CATALYST function for preparing the data
  # loading the beads data
  bc_key <- c(139, 141:156, 158:176)
  sce_bead <- assignPrelim(sce_bead, bc_key, verbose = FALSE)
  sce_bead <- applyCutoffs(estCutoffs(sce_bead))
  sce_bead <- computeSpillmat(sce_bead)
  counts_bead <- t(assay(sce_bead, "counts"))
  counts_bead <- floor(counts_bead)
  counts_bead <- as_tibble(counts_bead)
  
  channel_names <- rowData(sce_bead)[,"channel_name"]
  names(counts_bead) <- channel_names
  counts_bead <- mutate(counts_bead, barcode = sce_bead$bc_id)
  sm <- metadata(sce_bead)$spillover_matrix
  
  
  # loading the real cells data
  channel_names <- rowData(sce)[, "channel_name"]
  counts_real <- t(assay(sce, "counts"))
  counts_real <- floor(counts_real)
  colnames(counts_real) <- channel_names
  counts_real <- as_tibble(counts_real)
  
  
  # Starting the compensation procedure
  compensate <- function(counts_real, counts_bead, target_marker){
    
    spillover_markers <- names(which(sm[,target_marker] > 0))
    spillover_barcodes <- marker_to_code(spillover_markers)
    tb_bead <- counts_bead %>% 
      filter(barcode %in% spillover_barcodes) %>% 
      select(all_of(c(target_marker, "barcode"))) %>% 
      mutate(type = "beads")
    if(length(spillover_barcodes) > 1){
      # rename barcodes to marker names
      for(i in seq(spillover_barcodes)) {
        ids <- tb_bead$barcode == spillover_barcodes[i]
        tb_bead[ids, "barcode"] <- spillover_markers[i]
      }
    }
    
    tb_real <- counts_real %>% 
      select(all_of(target_marker)) %>%
      mutate(barcode = "none") %>%
      mutate(type = "real cells")
    # --------- mixture method ---------
    spillover_markers <- setdiff(spillover_markers, target_marker)
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
        filter(barcode == marker) %>% 
        pull(all_of(target_marker)) %>% 
        denoise(y_min = y_min, y_max = y_max) %>% 
        select(y, pmf)
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
      select(y, pmf)
    names(tb_real_pmf) <- c("y", target_marker)
    
    # join beads and real
    tb_pmf <- left_join(tb_beads_pmf, tb_real_pmf, by = "y")
    
    # --------- step 2: iterate ---------
    
    n_iter <- 10
    convergence <- matrix(nrow = n_iter, ncol = length(pi)+1)
    colnames(convergence) <- c("iteration", names(pi))
    
    for(i in 1:n_iter) {
      # remove counts without any signal
      total <- rowSums(select(tb_pmf, all_of(all_markers)))
      tb_pmf %<>% mutate(total = total)
      tb_pmf %<>% filter(total > 0)
      tb_pmf %<>% select(-total)
      
      # E-step
      
      # membership probabilities
      M <- tb_pmf %>% select(all_of(all_markers)) %>% as.matrix()
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
          filter(.data[[target_marker]] %in% ys) %>%
          pull(all_of(target_marker)) %>% 
          denoise(y_min = y_min, y_max = y_max) %>% 
          select(y, pmf)
        
      } else {
        
        # if no signal, then use uniform distribution  
        tb_real_pmf <- tb_beads_pmf %>% 
          select(y) %>% 
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
    M <- tb_pmf %>% select(all_of(all_markers)) %>% as.matrix()
    PI <- matrix(rep(pi, each = nrow(M)), ncol = ncol(M))
    post_M <- PI * M
    post_M <- post_M/rowSums(post_M)
    spill_prob <- 1-post_M[,target_marker]
    tb_spill_prob <- select(tb_pmf, y) %>% mutate(spill_prob = spill_prob)
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
    tb_compensate %<>% mutate(corrected = (1-spill)*.data[[target_marker]])
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
    res$target_marker <- target_marker
    class(res) <- "spillr"
    res
  }
  fit <- lapply(rownames(sm), function(marker) compensate(counts_real, counts_bead, marker))
  names(fit) <- rownames(sm)
  fit

}