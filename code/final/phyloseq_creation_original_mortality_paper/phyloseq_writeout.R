# write out count and phenotype data (handled in phyloseq objects in these scripts) to csvs
# see also functional_writeout.R

library(phyloseq)
library(microbiome)
library(dplyr)

source("code/final/phyloseq_creation/pseq_data_utilities.R")

pseq_phylum <- readRDS(phylum_level_phyloseq_path)
pseq_genus <- readRDS(genus_level_phyloseq_path)
pseq_species <- readRDS(species_level_phyloseq_path)


# the variables written out
variables_out <- c("Barcode", "BL_AGE", "MEN", "EAST", "BMI",
                   "CURR_SMOKE", "PREVAL_DIAB", "BL_USE_RX_L",
                   "BL_USE_RX_J01", "SYSTM", "BP_TREAT", "HFC_score", "UCOD",
                   "DEATH", "DEATH_AGEDIFF")

pheno_out <- meta(pseq_genus)[, variables_out]

# EAST, MEN are named factors - > recode to numeric (west=0, female=0)
pheno_out_num <- map_MEN_EAST_factors_to_num(pheno_out)
# ensure 0/1 variables are uniformly specified as factors
pheno_out_num$EAST <- as.factor(pheno_out_num$EAST)
pheno_out_num$MEN <- as.factor(pheno_out_num$MEN)
pheno_out_num$DEATH <- as.factor(pheno_out_num$DEATH)
pheno_out_num$CURR_SMOKE <- as.factor(pheno_out_num$CURR_SMOKE)
pheno_out_num$PREVAL_DIAB <- as.factor(pheno_out_num$PREVAL_DIAB)
pheno_out_num$BL_USE_RX_L <- as.factor(pheno_out_num$BL_USE_RX_L)
pheno_out_num$BL_USE_RX_J01 <- as.factor(pheno_out_num$BL_USE_RX_J01)
pheno_out_num$BP_TREAT <- as.factor(pheno_out_num$BP_TREAT)



outdir_pheno <- "output/fr02_mortality_datapublication/phenotype"

if (!dir.exists(outdir_pheno)) {
  dir.create(outdir_pheno, recursive=TRUE)
}

outdir_shogun <- "output/fr02_mortality_datapublication/shogun_tables"

if (!dir.exists(outdir_shogun)) {
  dir.create(outdir_shogun, recursive=TRUE)
}



#write.csv(pheno_out, file="input/data_work/phenotype_named.csv", row.names=TRUE)
write.csv(pheno_out_num, file=paste0(outdir_pheno, "/phenotype.csv"), row.names=TRUE)
write.csv(abundances(pseq_genus), file=paste0(outdir_shogun, "/table_genus.csv"), row.names=TRUE)
write.csv(abundances(pseq_species), file=paste0(outdir_shogun, "/table_species.csv"), row.names=TRUE)
write.csv(abundances(pseq_phylum), file=paste0(outdir_shogun, "/table_phylum.csv"), row.names=TRUE)
