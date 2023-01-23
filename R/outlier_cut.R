outlier_cut <- function(y, p ){
  quant_threshold <- quantile(y,probs= p)
  y[y<quant_threshold]
}
