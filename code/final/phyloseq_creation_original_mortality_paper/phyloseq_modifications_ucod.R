# replace K_TPKS variable (from registers) with aggregated variable UCOD
#  (underlying cause of death), with categories:
# C cancer
# I cardiovascular
# K gastrointestinal
# G neurological
# J respiratory
# O other

# also note that cause_of_death_enterobacteriaceae.R (as of writing this) assumes the old K_TPKS
# and extracts the equivalent contents with function get_causes_of_death

# I have added a new function get_causes_of_death_ucod that returns the UCOD variable created
# here; data is equivalent.

# Aaro Salosensaari May 2020

library(biomformat)
library(phyloseq)
library(data.table)
library(dplyr)
library(microbiome)
library(reshape2)

library(ggplot2)


add_ucod_to_pseq <- function(pseq) {
    # see also get_causes_of_death in R_functions.R
    causes_old <- meta(pseq)[, "K_TPKS"]
    causes <- gsub("[0-9]+", "", causes_old) # causes for all deaths

    causes_agg <- causes
    # aggregate all not CIKGJ into O, other
    for(i in 1:length(causes)) {
        if (!is.na(causes[i]) && !(causes[i] %in% c("C","I", "K", "G", "J"))) {
            causes_agg[i] <- "O"
        }
    }

    UCOD <- causes_agg
    sample_data(pseq) <- cbind(meta(pseq), UCOD)
    return(pseq)

}

# phyloseqs to be edited

phylum_level_phyloseq_noplasmid_hfc_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_phylum_ok_drop50k_2018-12-21_nop_hfc.rds"
genus_level_phyloseq_noplasmid_hfc_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_genus_all_drop50k_2018-12-21_nop_hfc.rds"
species_level_phyloseq_noplasmid_hfc_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_species_all_drop50k_2018-12-21_nop_hfc.rds"

pseq_phy <- readRDS(phylum_level_phyloseq_noplasmid_hfc_path)
pseq_phy <- add_ucod_to_pseq(pseq_phy)

pseq_gen <- readRDS(genus_level_phyloseq_noplasmid_hfc_path)
pseq_gen <- add_ucod_to_pseq(pseq_gen)

pseq_spe <- readRDS(species_level_phyloseq_noplasmid_hfc_path)
pseq_spe <- add_ucod_to_pseq(pseq_spe)

# new locations, no K_TPKS

phylum_level_phyloseq_noplasmid_hfc_ucod_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_phylum_ok_drop50k_2018-12-21_nop_hfc_ucod.rds"
genus_level_phyloseq_noplasmid_hfc_ucod_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_genus_all_drop50k_2018-12-21_nop_hfc_ucod.rds"
species_level_phyloseq_noplasmid_hfc_ucod_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_species_all_drop50k_2018-12-21_nop_hfc_ucod.rds"

saveRDS(pseq_phy, file=phylum_level_phyloseq_noplasmid_hfc_ucod_path)
saveRDS(pseq_gen, file=genus_level_phyloseq_noplasmid_hfc_ucod_path)
saveRDS(pseq_spe, file=species_level_phyloseq_noplasmid_hfc_ucod_path)