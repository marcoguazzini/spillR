extract_freq <- function(target, counts, threshold = 0.995, cut = FALSE){
  if(!cut){
x <- counts[,target]
tb_freq <- tibble(x) %>% group_by(x) %>% tally()
tb_freq$x %<>% as.integer()
tb_freq
}
  else{
    x <- counts[,target]
    x <- outlier_cut(x, p = threshold )
    tb_freq <- tibble(x) %>% group_by(x) %>% tally()
    tb_freq$x %<>% as.integer()
    tb_freq
  
  }
}
