library(biomformat)
library(phyloseq)
library(ggplot2)
library(vegan)
library(microbiome)

# Copypaste'd setup things here for now. not pretty, but works (easier to run with grun.py)

analysis_date <- paste0(Sys.Date())

pltdir <- paste0('output/figures/final/figure_ordination/', analysis_date, '/')
if (!dir.exists(pltdir)) {
  dir.create(pltdir, recursive=TRUE)
}

#datadir <- "output/raw_output/"

# Copypaste ends.


# Compute Bray-Curtis dissimilarity matrix, on species level taxonomic

phfs.pseq <- readRDS(species_level_phyloseq_path)


phfs.rel <- microbiome::transform(phfs.pseq, 'compositional')
otu <- abundances(phfs.rel)
bray.dist.m <- vegdist(t(otu), method="bray")

saveRDS(bray.dist.m, file=paste0(datadir, "phf_species_d50k_bray_dist_m",
                                 analysis_date, ".rds"))