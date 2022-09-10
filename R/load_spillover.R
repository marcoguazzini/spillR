load_spillover <- function(){
library(CATALYST)
library(SingleCellExperiment)
data(ss_exp, package = "CATALYST")
# specify mass channels stained for & debarcode
bc_ms <- c(139, 141:156, 158:176)
sce <- prepData(ss_exp)
sce <- assignPrelim(sce, bc_ms, verbose = FALSE)
sce <- applyCutoffs(estCutoffs(sce))
# compute & extract spillover matrix
sce <- computeSpillmat(sce)
sm <- metadata(sce)$spillover_matrix
# need to complete with zeros the spillover matrix
channel_names <- as.matrix(rowData(sce))[,"channel_name"]
smc <- matrix(data = 0, nrow = length(channel_names), ncol = length(channel_names))
rownames(smc) <- channel_names
colnames(smc) <- channel_names
for(i in rownames(sm)) {
  for(j in colnames(sm)) {
    smc[i,j] <- sm[i,j]
  }
}
diag(smc) <- 1.0
return(smc)
}
