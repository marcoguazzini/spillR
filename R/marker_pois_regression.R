marker_pois_regression <-
  function(target,
           barcode_target,
           barcode_emit,
           counts_spill,
          n_degree,
          smc) {
    
    tb_freq <- emit_extraction(target, barcode_target, barcode_emit, counts_spill, smc)
    # need to make sure to use raw = TRUE
    fit <- glm(
      n ~ poly(x, degree = n_degree, raw = TRUE),
      data = tb_freq,
      family = poisson(link = "log")
    )
    unname(fit[["coefficients"]])
  }
