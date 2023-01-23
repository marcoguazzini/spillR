outlier_cut <- function(y, prob ){
  quant_threshold <- quantile(y,probs = prob)
  y[y<quant_threshold]
}
