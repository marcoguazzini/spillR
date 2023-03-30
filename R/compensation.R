compensation <- function(target, fit_mix, tb_spills_coeff, y) {
  probs <- fit_mix@prior
  if (length(probs) < (length(tb_spills_coeff) + 1)) {
    zeros <- rep(0, ((length(
      tb_spills_coeff
    ) + 1) - length(probs)))
    probs <- c(probs, zeros)
  }
  marker_selection <- names(smc[, target])[smc[, target] > 0]
  tb_sel <- tibble(marker_name = marker_selection)
  tb_sel %<>% filter(marker_name != target)
  tb_sel$p_hat <- probs[-1]
  betas_target <-
    unname(fit_mix@components[["Comp.1"]][[1]]@parameters[["coef"]])
  tb_sel$betas <-
    sapply(tb_sel$marker_name, function(i)
      list(tb_spills_coeff[[i]]))
  tb_sel %<>% add_row(
    marker_name = target,
    p_hat = fit_mix@prior[1],
    betas = list(betas_target)
  )
  lambda_hat <- function(x, betas) {
    xp <- c(1, matrix(poly(x, degree = n_degree, raw = TRUE)))
    exp(sum(betas * xp))
  }
  
  calc_spill_prob <- function(x, tb_sel) {
    tb_hat <- tb_sel
    tb_hat %<>% mutate(x = x)
    tb_hat$lambda_hat <-
      sapply(1:nrow(tb_hat), function(i)
        lambda_hat(x, tb_hat$betas[[i]]))    ## bug the code returns Warning: longer object length is not a multiple of shorter object length but I am not getting why
    tb_hat %<>% mutate(pmf = dpois(x, lambda = lambda_hat))
    tb_hat %<>% mutate(spill_prob = p_hat * pmf)
    tb_hat %<>% mutate(spill_prob_norm = spill_prob / sum(tb_hat$spill_prob))
    
  }
  tb_prob <- lapply(unique(y), function(x)
    calc_spill_prob(x = x, tb_sel = tb_sel)) %>%
    bind_rows()
  
  tb <- tb_prob %>% filter(marker_name != target)
  tb <-
    data.frame(
      marker_name = tb$marker_name,
      x = tb$x,
      spill_prob_norm = tb$spill_prob_norm
    )
  tb_tot <-
    tb %>% group_by(x) %>% summarise(spill_prob = sum(spill_prob_norm)) %>%
    as.data.frame()
  n <- tibble(y) %>% group_by(y) %>% tally()
  tb_tot$n <- n$n
  tb_tot[is.na(tb_tot)] <- 0
  compl_prob <- function(p) {
    cp <- 1 - p
    if (cp < 0) {
      cp <- 0
    } else{
      cp
    }
  }
  tb_tot$c_prob <-
    sapply(tb_tot$spill_prob, function(p)
      compl_prob(p))
  expectation <- function(n, prob) {
    as.integer(n * prob)
  }
  
  tb_tot %<>% mutate(n_hat = expectation(n, c_prob))
  create_y <- function(y, n, n_hat) {
    if (n_hat > 0) {
      if (n - n_hat > 0) {
        comp_count <- c(rep(y, n_hat), rep(0, n - n_hat))
      }
      else{
        comp_count <- rep(y, n_hat)
      }
    } else{
      comp_count <- rep(0, n)
    }
    comp_count
  }
  comp <- unlist(sapply(1:nrow(tb_tot), function(i)
    create_y(tb_tot$x[i], tb_tot$n[i], tb_tot$n_hat[i])))
  sample(comp)
}
