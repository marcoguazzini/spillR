extract_freq <- function(target, counts, threshold , cut = FALSE){
  if(!cut){
x <- counts[,target]
tb_freq <- tibble(x) %>% group_by(x) %>% tally()
tb_freq$x %<>% as.integer()
tb_freq
}
  else{
    x <- counts[,target]
    x <- outlier_cut(x, threshold )
    tb_freq <- tibble(x) %>% group_by(x) %>% tally()
    tb_freq$x %<>% as.integer()
    tb_freq
  
  }
}
