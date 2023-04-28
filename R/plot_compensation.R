plot_compensation <- function(tb_compensate,target_marker){
  tfm <- function(x) asinh(x/5)
  colorscale = scale_fill_gradientn(
    colors = rev(brewer.pal(9, "YlGnBu")),
    values = c(0, exp(seq(-5, 0, length.out = 100)))
  )
    tb_compensate %>% 
      pivot_longer(c(all_of(target_marker), corrected), 
                   names_to = "status", values_to = target_marker) %>% 
      ggplot(aes(tfm(.data[[target_marker]]), fill = status)) + 
      geom_histogram(position = "identity", alpha = 0.3, bins = 50) + 
      facet_wrap(~status, ncol = 1)
  }
