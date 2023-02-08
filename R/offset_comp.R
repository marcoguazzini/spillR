offset_comp <- function(tb_freq,betas,target, barcode_target, barcode_emit, counts_spill,smc){
  tb <- emit_extraction(target, barcode_target, barcode_emit,smc)
  offset <- function(x,betas,tb){
    if(x <= max(range(tb$x)) && x >= min(range(tb$x))){
      xp <- c(1, matrix(poly(x, degree = n_degree, raw = TRUE)))
      sum(betas * xp)
    }
    else{
      -10^7
    }
  }
  unlist(lapply(tb_freq$x, function(xi) offset(xi,betas, tb)))
}
