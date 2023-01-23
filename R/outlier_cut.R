outlier_cut <- function(y, prob ){
  quant_threshold <- quantile(y,prob)
  y[y<quant_threshold]
}
