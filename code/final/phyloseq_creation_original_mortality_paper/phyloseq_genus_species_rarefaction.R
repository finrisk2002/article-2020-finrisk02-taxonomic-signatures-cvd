# rarefy for checks
# rarefaction done by sampling without replacement for least worse rarefaction

library(biomformat)
library(phyloseq)
library(data.table)
library(dplyr)
library(microbiome)
library(reshape2)


genus_level_phyloseq_noplasmid_hfc_ucod_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_genus_all_drop50k_2018-12-21_nop_hfc_ucod.rds"

genus_level_phyloseq_noplasmid_hfc_ucod_rar_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_genus_all_drop50k_2018-12-21_nop_hfc_ucod_rar10pctl.rds"


species_level_phyloseq_noplasmid_hfc_ucod_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_species_all_drop50k_2018-12-21_nop_hfc_ucod.rds"

species_level_phyloseq_noplasmid_hfc_ucod_rar_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_species_all_drop50k_2018-12-21_nop_hfc_ucod_rar10pctl.rds"


pseq_g <- readRDS(genus_level_phyloseq_noplasmid_hfc_ucod_path)
pseq_s <- readRDS(species_level_phyloseq_noplasmid_hfc_ucod_path)


gs <- quantile(sample_sums(pseq_g), probs=c(0.10))
pseq_g_rar <- rarefy_even_depth(pseq_g, sample.size = gs[[1]], replace=FALSE, rngseed=1234)



ss <- quantile(sample_sums(pseq_s), probs=c(0.10))
pseq_s_rar <- rarefy_even_depth(pseq_s, sample.size = ss[[1]], replace=FALSE, rngseed=1234)


saveRDS(pseq_g_rar, file=genus_level_phyloseq_noplasmid_hfc_ucod_rar_path)
tmsaveRDS(pseq_s_rar, file=species_level_phyloseq_noplasmid_hfc_ucod_rar_path)
