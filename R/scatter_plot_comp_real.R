scatter_plot_comp_real <-  function(chs, counts_nnls,counts, y_comp_mix){
  colorscale = scale_fill_gradientn(
    colors = rev(brewer.pal(9, "YlGnBu")),
    values = c(0, exp(seq(-5, 0, length.out = 100)))
  )
  tfm <- function(x) asinh(x/5)
  # spillover without compensation
  tb_counts_uncomp <- counts %>% as_tibble() %>% dplyr::select(all_of(chs))
  tb_counts_uncomp %<>% mutate(cell = 1:nrow(tb_counts_uncomp))
  tb_counts_uncomp %<>% mutate(comp = "with spillover")
  # compensate using NNLS as in CATALYST::compCytof
  tb_counts_nnls <- counts_nnls %>% as_tibble() %>% dplyr::select(all_of(chs))
  tb_counts_nnls %<>% mutate(cell = 1:nrow(tb_counts_nnls))
  tb_counts_nnls %<>% mutate(comp = "non-negative least squares")
  
  # generate with our method spillR::compensate
  tb_counts_polyspillr <- y_comp_mix %>% as_tibble() %>% dplyr::select(all_of(chs))
  tb_counts_polyspillr %<>% mutate(cell = 1:nrow(tb_counts_polyspillr))
  tb_counts_polyspillr %<>% mutate(comp = "our method spillR")
  # merge and plot
  tb_counts_combined <- bind_rows(
                                  tb_counts_uncomp, 
                                  tb_counts_nnls, 
                                  tb_counts_polyspillr)
  tb_counts_combined$comp %<>% factor(
    levels = c( 
               "with spillover", 
               "non-negative least squares",
               "our method spillR")
  )
  
  p1 <- ggplot(tb_counts_combined, aes_string(x = chs[1], y = chs[2])) +
    geom_hex(bins = 64) +
    coord_fixed() +
    colorscale + 
    facet_wrap(~comp, nrow = 1) 
  p2 <-  ggplot(tb_counts_combined %>% mutate_at(chs, tfm), 
                aes_string(x = chs[1], y = chs[2])) +
    geom_hex(bins = 64) +
    coord_fixed() +
    colorscale + 
    facet_wrap(~comp, nrow = 1)
  
  list(p1,p2)
}