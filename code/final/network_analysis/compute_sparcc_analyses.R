# compute spieceasi analyses
source("code/final/R_functions/isinstalled.R")
source("code/final/R_functions/taxa_names_functions.R")

is.installed(c("biomformat", "microbiome", "phyloseq", "data.table",
			   "RColorBrewer", "dplyr", "reshape2", "SpiecEasi"))


phfg.pseq.d50k <- readRDS(genus_level_phyloseq_path)

pseq.d50k.rel.core <- core(microbiome::transform(phfg.pseq.d50k, "compositional"),
						   detection=0.1/100, prevalence=1/100)

pseq.d50k.core <- prune_taxa(taxa_names(pseq.d50k.rel.core), phfg.pseq.d50k)

sparcc_res <- sparcc(t(abundances(pseq.d50k.core)))

saveRDS(sparcc_res, paste0(datadir,"/sparcc_pseq_2018-12-21_drop50k_core_genus_glasso_defaults.RDs"))