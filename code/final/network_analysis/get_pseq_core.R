
get_pseq_genus_core <- function() {
	if(!exists("genus_level_phyloseq_path")) {
		cat("load main_setup.R?\n")
	}
	pseq <- readRDS(genus_level_phyloseq_path)

	# pretty names
	taxa_names(pseq) <- gsub("^g_", "", taxa_names(pseq))

	pseq.rel.core <- core(microbiome::transform(pseq, "compositional"), detection=0.1/100, prevalence=1/100)

	pseq.core <- prune_taxa(taxa_names(pseq.rel.core), pseq)
}