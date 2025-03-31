# add / make new / check existing covariates for CVD analyses


library(biomformat)
library(phyloseq)
library(data.table)
library(dplyr)
library(microbiome)
library(reshape2)

# original paths

# NB. needs to point to directory where main input data files are located
input_datadir <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/"
#input_datadir <- "input/data_main/"

# Paths for the phyloseq objects, saved in .RDS format
phylum_level_phyloseq_noplasmid_hfc_ucod_path <- paste0(input_datadir,
    "phfinrisk_phylum_ok_drop50k_2018-12-21_nop_hfc_ucod.rds")
genus_level_phyloseq_noplasmid_hfc_ucod_path <- paste0(input_datadir,
    "phfinrisk_genus_all_drop50k_2018-12-21_nop_hfc_ucod.rds")
species_level_phyloseq_noplasmid_hfc_ucod_path <- paste0(input_datadir,
    "phfinrisk_species_all_drop50k_2018-12-21_nop_hfc_ucod.rds")


# covariates, check that are included

# BL_AGE, MEN, non-HDL COL (= KOL - HDL), CURR_SMOKE, SYSTM, PREVAL_DIAB, BMI


# other notes: in survival analysis
# - remove cases with prevalent CVD (PREVAL_CVD)
# - analyse prevalent CVD (INCIDENT_CVD, CVD_AGEDIFF)

# updated paths

phylum_level_phyloseq_path <-  paste0(input_datadir,
    "phfinrisk_phylum_ok_drop50k_2018-12-21_nop_hfc_ucod_cvd.rds")
genus_level_phyloseq_path <- paste0(input_datadir,
    "phfinrisk_genus_all_drop50k_2018-12-21_nop_hfc_ucod_cvd.rds")
species_level_phyloseq_path <- paste0(input_datadir,
    "phfinrisk_species_all_drop50k_2018-12-21_nop_hfc_ucod_cvd.rds")
