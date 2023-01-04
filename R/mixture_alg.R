mixture_alg <- function(target, x, barcode_target, marker_selection, bar, thr = 0.995, cut = FALSE){
  tb_spills_coeff <- lapply(bar, function(b) marker_pois_regression(target = target, barcode_target = barcode_target, barcode_emit=b, counts_spill = counts_spill))
  tb_freq <- extract_freq(target, counts,threshold = thr,  cut = cut)
  fit_mix <- fit_mixture(tb_freq, tb_spills_coeff, target,barcode_target,bar,counts_spill)
  compensation(target, fit_mix, tb_spills_coeff,x)
}