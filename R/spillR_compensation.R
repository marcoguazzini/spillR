spillR_compensation <- function(markers, counts_real, counts_bead, n_iter,sm){
  
  EM_mixture <- function(target_marker, counts_real, counts_bead, n_iter,sm){
    marker_to_code <- function(marker) as.integer(substr(marker, 3, 5))
    target_barcode <- marker_to_code(target_marker)
    tfm <- function(x) asinh(x/5)
    colorscale = scale_fill_gradientn(
      colors = rev(brewer.pal(9, "YlGnBu")),
      values = c(0, exp(seq(-5, 0, length.out = 100)))
    )
    
    denoise <- function(y, degree, y_max = max(y) ) {
      
      # frequency table
      tb_obsv <- tibble(y) %>% group_by(y) %>% tally()
      y_min_obsv <- min(tb_obsv$y)
      y_max_obsv <- max(tb_obsv$y)
      tb_pred <- tibble(y = 0:y_max)
      tb_pred %<>% left_join(tb_obsv, by = "y")
      
      # padding with zero outside of data support
      tb_pred %<>% mutate(n = ifelse(y > y_max_obsv | y < y_min_obsv, 0, n))
      
      # option 4
      tb_pred$n[is.na(tb_pred$n)] <- 0
      tb_pred %<>% mutate(lambda = n)
      
      # normalize
      tb_pred %>% mutate(pmf = lambda/sum(tb_pred$lambda))
      
    }
    
    # creating tb_real
    tb_real <- counts_real %>% 
      select(all_of(target_marker)) %>%
      mutate(barcode = "none") %>%
      mutate(type = "real cells")
    
    # what markers spillover into target marker?
    spillover_markers <- names(which(sm[,target_marker] > 0))
    spillover_barcodes <- marker_to_code(spillover_markers)
    tb_bead <- counts_bead %>% 
      filter(barcode %in% spillover_barcodes) %>% 
      select(all_of(c(target_marker, "barcode"))) %>% 
      mutate(type = "beads")
    
    # getting y_max
    y_max <- tb_real %>% pull(all_of(target_marker)) %>% max()
    
    # support for target marker
    tb_beads_pmf <- tibble(y = 0:y_max)
    
    # collect pmf from beads
    for(marker in setdiff(spillover_markers, target_marker)) {
      
      tb <- tb_bead %>% 
        filter(barcode == marker_to_code(marker)) %>% 
        pull(all_of(target_marker)) %>% 
        denoise(degree = degree_bead, y_max = y_max) %>% 
        select(y, pmf)
      names(tb) <- c("y", marker)
      tb_beads_pmf %<>% left_join(tb, by = "y")
      
    }
    
    
    # prior probability
    pi <- rep(1, length(spillover_markers))
    pi <- pi/length(pi)
    names(pi) <- spillover_markers
    
    # add pmf from real cells
    tb_real_pmf <- tb_real %>% 
      pull(all_of(target_marker)) %>% 
      denoise(degree = degree_real, y_max = y_max) %>%
      select(y, pmf)
    names(tb_real_pmf) <- c("y", target_marker)
    
    # join beads and real
    tb_pmf <- left_join(tb_beads_pmf, tb_real_pmf, by = "y")
    
    # ALGORITHM
    convergence <- matrix(nrow = n_iter, ncol = length(pi)+1)
    colnames(convergence) <- c("iteration", names(pi))
    
    for(i in 1:n_iter) {
      # remove counts without any signal
      total <- rowSums(select(tb_pmf, all_of(spillover_markers)))
      tb_pmf %<>% mutate(total = total)
      tb_pmf %<>% filter(total > 0)
      tb_pmf %<>% select(-total)
      
      # E-step
      
      # membership probabilities
      M <- tb_pmf %>% select(all_of(spillover_markers)) %>% as.matrix()
      PI <- matrix(rep(pi, each = nrow(M)), ncol = ncol(M))
      post_M <- PI * M
      post_M <- post_M/rowSums(post_M)
      
      # M-step
      
      # update prior probability
      pi <- colSums(post_M)/nrow(post_M)
      
      # stochastic EM: 
      # assigns each observation to a class with the highest posterior probability
      #class <- apply(M, 1, function(Mi) sample(colnames(M), size = 1, prob = Mi))
      # categorical EM:
      # assigns each observation randomly based on posterior probabilities
      class <- spillover_markers[apply(M, 1, which.max)]
      
      # update pmf from real cells
      if(sum(class == target_marker) > 0) {
        
        ys <- tb_pmf[class == target_marker, ] %>% pull(y)
        tb_real_pmf <- tb_real %>% 
          filter(.data[[target_marker]] %in% ys) %>%
          pull(all_of(target_marker)) %>% 
          denoise(degree = degree_real, y_max = y_max) %>%
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
    
    # calculate posterior spillover probability for each cell
    tb_pmf <- round(tb_pmf, digits = 6)
    M <- tb_pmf %>% select(all_of(spillover_markers)) %>% as.matrix()
    PI <- matrix(rep(pi, each = nrow(M)), ncol = ncol(M))
    post_M <- PI * M
    post_M <- post_M/rowSums(post_M)
    spill_prob <- 1-post_M[,target_marker]
    tb_spill_prob <- select(tb_pmf, y) %>% mutate(spill_prob = spill_prob)
    tb_spill_prob %<>% mutate(spill_prob = if_else(is.na(spill_prob), 0, spill_prob))
    names(tb_spill_prob) <- c(target_marker, "spill_prob")
    
    
    # compensate
    tb_compensate <- data.frame(count= unlist(counts_real[,target_marker]))
    colnames(tb_compensate)<- target_marker
    tb_compensate %<>% left_join(tb_spill_prob, by = target_marker)
    tb_compensate %<>% mutate(
      spill = rbinom(n = nrow(tb_compensate), 
                     size = 1, 
                     prob = tb_compensate$spill_prob)
    )
    tb_compensate %<>% mutate(corrected = (1-spill)*.data[[target_marker]])
    tb_spill_prob %<>% mutate(y_tfm = tfm(.data[[target_marker]]))
    fit <- glm(spill_prob ~ y_tfm, family = binomial, data = tb_spill_prob)
    inverse_logit <- function(fit, x) {
      hat <- coef(fit)[1] + coef(fit)[2]*x
      1/(1+exp(-hat))
    }
    tb_spill_prob %<>% mutate(spill_prob_smooth = inverse_logit(fit, y_tfm))
    list(tb_compensate, tb_spill_prob)
  }
  fit<-
    lapply(markers, function(marker)
      EM_mixture(
        marker,
        counts_real,
        counts_bead,
        n_iter = 20,
        sm
      )
    )
  fit
}