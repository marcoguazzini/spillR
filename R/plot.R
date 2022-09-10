plot <- function(counts, comp_counts, chan_name){
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(magrittr)
  library(cowplot)
  library(flowCore)
  library(ggplot2)
  
  ch1 <- chan_name[1]
  ch2 <- chan_name[2]
  
  colorscale = scale_fill_gradientn(
    colors = rev(brewer.pal(9, "YlGnBu")),
    values = c(0, exp(seq(-5, 0, length.out = 100)))
  )
  tfm <- function(x) asinh(x/5)
  tb_counts <- as_tibble(counts)
  
  # original
  tb_counts_uncomp <- tb_counts %>% select(ch1, ch2)
  tb_counts_uncomp %<>% mutate(cell = 1:nrow(tb_counts))
  tb_counts_uncomp %<>% mutate(comp = "no compensation")
  
  # compensate
  tb_counts_comp <- as_tibble(comp_counts)
  tb_counts_comp <- tb_counts_comp %>% select(ch1, ch2)
  #tb_counts_comp %<>% mutate(cells1 = Er167Di-10)
  #tb_counts_comp %<>% mutate(cells2 = Er168Di-10)
  tb_counts_comp %<>% mutate(cell = 1:nrow(tb_counts))
  tb_counts_comp %<>% mutate(comp = "our method spillR")
  
  tb_counts_combined <- bind_rows(tb_counts_uncomp, tb_counts_comp)
  tb_counts_combined <- bind_rows(tb_counts_uncomp, tb_counts_comp)
  tb_counts_combined_transformed <- tb_counts_combined %>% mutate( tfm(tb_counts_combined[ch1]))
  tb_counts_combined_transformed  %<>% mutate( tfm(tb_counts_combined_transformed[ch2]))
  tb_counts_combined_transformed %<>% select(ch1,ch2,cell,comp)
  
plot_grid(
 ggplot(tb_counts_combined, aes_string(x = ch1, y =ch2)) +
    geom_hex(bins = 64) +
    coord_fixed() +
    colorscale +
    facet_wrap(~comp),
  
 ggplot(tb_counts_combined_transformed, aes_string(x = ch1, y = ch2)) +
    geom_hex(bins = 64) +
    coord_fixed() +
    colorscale +
    facet_wrap(~comp))
 
}
