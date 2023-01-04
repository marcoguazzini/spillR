plot_compensation <- function(y, y_comp){
  tb_comp <- data.frame(y = y_comp)
  tb_comp %<>% mutate( comp = "compensated")
  tb_uncomp <- data.frame(y = y)
  tb_uncomp %<>% mutate(comp = "uncompensated")
  tb_counts_combined <- bind_rows(tb_uncomp, tb_comp)
  
  
  ggplot(tb_counts_combined, aes(x=y,fill=comp)) +
    geom_histogram(alpha = 0.3, position = 'identity', bins = 50)+
    theme(legend.position="top")+
    scale_fill_manual(values= c("blue", "green"))}