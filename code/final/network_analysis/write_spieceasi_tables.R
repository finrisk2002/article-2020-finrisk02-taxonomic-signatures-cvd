library(SpiecEasi)

# also write SpiecEasi "correlations" as CSV

if (!exists("pseq.core")) {
	source("code/final/network_analysis/get_pseq_core.R")
	pseq.core <- get_pseq_genus_core()
}


spiec.res <- readRDS(paste0(datadir,"/spieceasi_pseq_2018-12-21_drop50k_core_genus_glasso_nl30_th001.RDs"))
# adjacency / correlation matrix
spiec.res.cor  <- cov2cor(as.matrix(getOptCov(spiec.res)))
colnames(spiec.res.cor) <- rownames(spiec.res.cor) <- colnames(t(abundances(pseq.core)))

write.csv(spiec.res.cor, paste0(outdir, "table_shogun_spieceasi_correlation.csv"))