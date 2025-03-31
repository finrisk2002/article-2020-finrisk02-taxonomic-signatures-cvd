library(microbiome)
library(vegan)

source("code/final/create_figure1_map_and_pcoa/setup.R")

# Load the FINRISK data set
datadir <- "input/data_work/"
pcoa.file <- paste0(datadir, "phf_species_d50k_bray_pseq_pcoa2020-02-13.rds")

# SHOGUN
genus.file   <- genus_level_phyloseq_path
species.file   <- species_level_phyloseq_path

phfinrisk_genus <- readRDS(genus.file)
phfinrisk_species <- readRDS(species.file)
phfs.rel.ord.pcoa <- readRDS(pcoa.file)



