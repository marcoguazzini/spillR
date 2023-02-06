mixture_alg <- function(counts, smc, n_degree, t = 0.995, c = FALSE){
channel_names <- colnames(counts) 
targets <- unname(channel_names[-(1:4)])
barcode_targets <- unname(sapply(targets , function(t) substr(t, 3,5)))
marker_selections <- sapply(targets, function(t) list(names(smc[, t])[smc[, t] > 0]))
extr_bars <- function(markers, barcode_target){
  p <- sapply(markers, function(m) substr(m, 3, 5))
  p[p != barcode_target]
}
bar <- sapply(1:length(barcode_targets) , function(i) list(extr_bars(marker_selections[[i]], barcode_targets[i] )))

mixture <- function(target, x, barcode_target, marker_selection, bar, thr = t, cut = c, n = n_degree){
  tb_spills_coeff <- lapply(bar, function(b) marker_pois_regression(target = target, barcode_target = barcode_target, barcode_emit=b, n_degree = n, smc = smc))
  tb_freq <- extract_freq(target, counts, threshold = t,  cut = c)
  fit_mix <- fit_mixture(tb_freq, tb_spills_coeff, target,barcode_target,bar,counts_spill, n, smc)
  compensation(target, fit_mix, tb_spills_coeff,x)
}
counts_poly_mix <-data.frame(lapply(1:length(barcode_targets), function(i) mixture(targets[i], counts[,targets[i]],barcode_targets[i],marker_selections[[i]], bar[[i]])))
counts_poly_mix <- data.frame(counts[,1:4],counts_poly_mix)
colnames(counts_poly_mix ) <- channel_names
counts_poly_mix 
}
