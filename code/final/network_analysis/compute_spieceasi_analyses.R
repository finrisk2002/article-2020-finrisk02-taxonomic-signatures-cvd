# compute spieceasi analyses
source("code/final/R_functions/isinstalled.R")
source("code/final/R_functions/taxa_names_functions.R")

is.installed(c("biomformat", "microbiome", "phyloseq", "data.table",
			   "RColorBrewer", "dplyr", "reshape2", "SpiecEasi"))

source("code/final/network_analysis/pseq_spieceasi.R")

phfg.pseq.d50k <- readRDS(genus_level_phyloseq_path)

pseq.d50k.rel.core <- core(microbiome::transform(phfg.pseq.d50k, "compositional"),
						   detection=0.1/100, prevalence=1/100)

pseq.d50k.core <- prune_taxa(taxa_names(pseq.d50k.rel.core), phfg.pseq.d50k)

spieceasi_th001 <- pseq_spieceasi(pseq.d50k.core, lambda.min.ratio=1e-2, nlambda=30,
							  	  rep.num=50, seed=42, thresh=0.01)

# stability threshodl 0.01 used
saveRDS(spieceasi_th001, paste0(datadir,"/spieceasi_pseq_2018-12-21_drop50k_core_genus_glasso_nl30_th001.RDs"))

# stability threshold 0.05 not used
spieceasi_th005 <- pseq_spieceasi(pseq.d50k.core, lambda.min.ratio=1e-2, nlambda=30,
							  	  rep.num=50, seed=42, thresh=0.05)
saveRDS(spieceasi_th005, paste0(datadir,"/spieceasi_pseq_2018-12-21_drop50k_core_genus_glasso_nl30_th005.RDs"))

# east west replication at th001

pseq.d50k.core.east <- prune_samples(sample_names(pseq.d50k.core)[meta(pseq.d50k.core)[,"EAST"] == "EAST"],
                                                 pseq.d50k.core)
pseq.d50k.core.west <- prune_samples(sample_names(pseq.d50k.core)[meta(pseq.d50k.core)[,"EAST"] == "WEST"],
                                                 pseq.d50k.core)
spieceasi_th001_east <- pseq_spieceasi(pseq.d50k.core.east, lambda.min.ratio=1e-2, nlambda=30,
							  	  rep.num=50, seed=42, thresh=0.01)
saveRDS(spieceasi_th001_east,
	paste0(datadir,"/spieceasi_pseq_2018-12-21_drop50k_core_genus_glasso_nl30_th001_east.RDs"))

spieceasi_th001_west <- pseq_spieceasi(pseq.d50k.core.west, lambda.min.ratio=1e-2, nlambda=30,
							  	  rep.num=50, seed=42, thresh=0.01)
saveRDS(spieceasi_th001_west,
	paste0(datadir,"/spieceasi_pseq_2018-12-21_drop50k_core_genus_glasso_nl30_th001_west.RDs"))