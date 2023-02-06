emit_extraction <-
  function(target,
           barcode_target,
           barcode_emit,
          smc) {
    marker_selection <- names(smc[, target])[smc[, target] > 0]
    emit_spill <-
      sapply(marker_selection, function(marker)
        substr(marker, 3, 5))
 counts_spill <- load_counts_spill()
    y <- counts_spill %>%
      select_at(c(barcode_target, "barcode")) %>%
      filter(barcode %in% emit_spill) %>% filter(barcode != barcode_target)
    y %<>% filter(barcode == barcode_emit)
    colnames(y)[1] <- "x"
    x <- data.frame(x = unname(y$x))
    tb_freq <- tibble(x) %>% group_by(x) %>% tally()
    tb_freq$x %<>% as.integer()
    tb_freq
  }
