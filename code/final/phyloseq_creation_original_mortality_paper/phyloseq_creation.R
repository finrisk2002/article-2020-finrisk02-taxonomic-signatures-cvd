# this script contains notes about creation of the internal phyloseqs from Institute data sources (Aki Havulinna)
# intended to work on FIMM Atlas server only

library(biomformat)
library(phyloseq)
library(data.table)
library(dplyr)
library(microbiome)


source("code/final/phyloseq_creation/pseq_data_utilities.R")


## SHOGUN

# bioms from KnightLab
taxa_counts_bioms_dir <- "/csc/fr_metagenome/microbiome2018-02-26/taxonomy/"
genus_shogun_counts_file <- paste0(taxa_counts_bioms_dir, "shogun/combined_redist.genus.biom")
species_shogun_counts_file <- paste0(taxa_counts_bioms_dir, "shogun/combined_redist.species.biom")
phylum_shogun_counts_file <- paste0(taxa_counts_bioms_dir, "shogun/combined_redist.phylum.biom")

# paths to shogun phyloseqs that will be created
# date describes the script, not data.

phylum_level_phyloseq_orig_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_shogun_phylum_all_drop50k_2020-03-30.rds"
genus_level_phyloseq_orig_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_shogun_genus_all_drop50k_2020-03-30.rds"
species_level_phyloseq_orig_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_shogun_species_all_drop50k_2020-03-30.rds"

pheno_file <- "/csc/fr_metagenome/pheno/2015_60_Salomaa_Jain_dataFR02_2018-08-15.RData" 
pheno_data <- read_pheno_data2(pheno_file)

# covariate descriptions, csv
desc.file <- "input/covariate_descriptions.csv"
phfinrisk_metadatadesc <- read_metadata_descriptions(desc.file)

# transform categorical variables to factors and other post-processing
pheno_data_clean <- clean_pheno_data(pheno_data, phfinrisk_metadatadesc)
# create shogun Phyloseq-object from biom files and cleaned data
phfinrisk_genus_all   <- process_shogun_biom(genus_shogun_counts_file, pheno_data_clean)
# repeat for phylum, species
phfinrisk_phylum_all <- process_shogun_biom(phylum_shogun_counts_file, pheno_data_clean)
phfinrisk_species_all <- process_shogun_biom(species_shogun_counts_file, pheno_data_clean)

# remove samples with < 50 k counts
phfinrisk_genus_d50k <- drop_small_rc_pseq(phfinrisk_genus_all, countthr=50000)
phfinrisk_species_d50k <- drop_small_rc_pseq(phfinrisk_species_all, countthr=50000)
phfinrisk_phylum_d50k <- drop_small_rc_pseq(phfinrisk_phylum_all, countthr=50000)

saveRDS(phfinrisk_genus_d50k, file = genus_level_phyloseq_orig_path)
saveRDS(phfinrisk_phylum_d50k, file = phylum_level_phyloseq_orig_path)
saveRDS(phfinrisk_species_d50k, file = species_level_phyloseq_orig_path)


## Centrifuge; similar, but update paths.

# genus_centr_counts_file <- paste0(taxo_dataraw_dir, "centrifuge/combined_profile.genus.biom")
# species_centr_counts_file <- paste0(taxo_dataraw_dir, "centrifuge/combined_profile.species.biom")
# phylum_centr_counts_file <- paste0(taxo_dataraw_dir, "centrifuge/combined_profile.phylum.biom")



# phfinrisk_centr_genus_all   <- process_centrifuge_biom(genus_centr_counts_file, pheno_data_clean)
# phfinrisk_centr_species_all   <- process_centrifuge_biom(species_centr_counts_file, pheno_data_clean)
# phfinrisk_centr_phylum_all   <- process_centrifuge_biom(phylum_centr_counts_file, pheno_data_clean)


# ### Centrifuge, supplemental
# # very slow way to extract prettier taxa names

# taxa_counts_reports_dir <- taxa_counts_bioms_dir
# allsubdirs <- list.dirs(taxa_counts_reports_dir)

# check_if_dir_contains_profile <- function(dirpath) {
#   # potential stoolid is that directory name
#   stoolid <- tail(strsplit(dirpath,'/')[[1]], n=1)
#   if (dir.exists(file.path(dirpath, "centrifuge"))) {
#     reportfile <- file.path(dirpath, "centrifuge", paste0(stoolid, ".report.txt"))
#     if (file.exists(reportfile)) {
#       return(c(dirpath, stoolid))
#     }
#   }
#   return(NULL)
# }

# reportdir_id_pairs <- lapply(allsubdirs, check_if_dir_contains_profile)

# reportdirs <- sapply(reportdir_id_pairs[!sapply(reportdir_id_pairs, is.null)],
#                           function (pair) {pair[[1]]})
# reportnames <- sapply(reportdir_id_pairs[!sapply(reportdir_id_pairs, is.null)],
#                           function (pair) {pair[[2]]})
# reportdirs_named <- reportdirs
# names(reportdirs_named) <- reportnames

# # create Phyloseq-object from biom files and cleaned data
# phfinrisk_centr_genus_all2   <- process_centrifuge_biom(genus_centr_counts_file, pheno_data_clean, reportdirs_named)
# phfinrisk_centr_species_all2   <- process_centrifuge_biom(species_centr_counts_file, pheno_data_clean, reportdirs_named)
# phfinrisk_centr_phylum_all2   <- process_centrifuge_biom(phylum_centr_counts_file, pheno_data_clean, reportdirs_named)