extract_freq <- function(target, counts, threshold , cut = FALSE){
  outlier_cut <- function(y, prob){
  quant_threshold <- quantile(y,probs = prob)
  y[y<quant_threshold]
}
  if(!cut){
x <- counts[,target]
tb_freq <- tibble(x) %>% group_by(x) %>% tally()
tb_freq$x %<>% as.integer()
tb_freq
}
  else{
    x <- counts[,target]
    x <- outlier_cut(x, prob =threshold )
    tb_freq <- tibble(x) %>% group_by(x) %>% tally()
    tb_freq$x %<>% as.integer()
    tb_freq
  
  }
}
