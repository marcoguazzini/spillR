#' Generate compensated counts dataset
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import RColorBrewer
#' @import magrittr
#' @import cowplot
#' @import flowCore
#' @import ggplot2
#' @export
#'
#' @param counts Matrix with counts, on the rows the cells, on the column the markers
#'
#' @param comp_counts Matrix with compensated counts, on the rows the cells, on the column the markers
#'
#' @param chan_names  Desired channels we want to compare
#'   
#' @return plot
#'
#' @examples
#' set.seed(23)
#' sm <- spillR::load_spillover()
#' counts <- spillR::prepare_data(sm, shape, rate, n_cells)
#' compensated_counts <- generate_data(counts,spillover_matrix)
#' spillR::plot(counts, comp_counts, chan_name)

plot <- function(counts, comp_counts, chan_name){

  
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
  ggtitle("Untransformed data")+
    coord_fixed() +
    colorscale +
    facet_wrap(~comp),
  
 ggplot(tb_counts_combined_transformed, aes_string(x = ch1, y = ch2)) +
    geom_hex(bins = 64) +
   ggtitle("Transformed data with asinh(x/5) transformation")+
    coord_fixed() +
    colorscale +
    facet_wrap(~comp))
}
