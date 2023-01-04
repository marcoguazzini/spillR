outlier_cut <- function(y, p = 0.995){
  quant_threshold <- quantile(y, probs = p)
  y[y<quant_threshold]
}
