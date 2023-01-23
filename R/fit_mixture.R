fit_mixture <-
  function(tb_freq,
           tb_spills_coeff,
           target,
           barcode_target,
           bar,
           counts_spill,
          n_degree,
           smc) {
    offset <-
      lapply(1:length(bar), function(i)
        offset_comp(
          tb_freq,
          tb_spills_coeff[[i]],
          target,
          barcode_target,
          bar[i],
          counts_spill,
          smc
        ))
    model <- list(FLXglm(n ~ poly(x, degree = n_degree, raw = TRUE), family = "poisson"))
    model <- append(model,
                    lapply(1:length(bar), function(i)
                      FLXglm(n ~ -1 , offset = offset[[i]], family = "poisson")))
    
    
    fit_mix <- stepFlexmix(~ x,
                           data = tb_freq,
                           k = length(bar) + 1,
                           model = model)
    return(fit_mix)
  }
