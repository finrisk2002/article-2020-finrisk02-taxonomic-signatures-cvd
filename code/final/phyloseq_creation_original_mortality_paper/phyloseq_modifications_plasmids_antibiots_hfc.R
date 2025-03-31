# phyloseq operations
# * remove plasmids
# * add HFC_score
# * write pseqs with and without rxj01 (antibiotics use) variable


library(biomformat)
library(phyloseq)
library(data.table)
library(dplyr)
library(microbiome)
library(reshape2)

library(ggplot2)

# remove Plasmids from Shogun

phylum_level_phyloseq_orig_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_phylum_ok_drop50k_2018-12-21.RDs"
genus_level_phyloseq_orig_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_genus_all_drop50k_2018-12-21.RDs"
species_level_phyloseq_orig_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_species_all_drop50k_2018-12-21.RDs"

pseq_phylum <- readRDS(phylum_level_phyloseq_orig_path)
pseq_genus <- readRDS(genus_level_phyloseq_orig_path)
pseq_species <- readRDS(species_level_phyloseq_orig_path)

phylum_plasmids <- grepl("Plasmid", taxa_names(pseq_phylum))
genus_plasmids <- grepl("Plasmid", taxa_names(pseq_genus))
species_plasmids <- grepl("Plasmid", taxa_names(pseq_species))

pseq_phylum_noplasmid <- prune_taxa(!phylum_plasmids, pseq_phylum)
pseq_genus_noplasmid <- prune_taxa(!genus_plasmids, pseq_genus)
pseq_species_noplasmid <- prune_taxa(!species_plasmids, pseq_species)

# investigate


pseq_phylum_plasmid <- prune_taxa(phylum_plasmids, pseq_phylum)
pseq_genus_plasmid <- prune_taxa(genus_plasmids, pseq_genus)
pseq_species_plasmid <- prune_taxa(species_plasmids, pseq_species)

species_df <- data.frame(Sample=sample_names(pseq_species),
                         Plasmid=sample_sums(pseq_species_plasmid),
                         Other=sample_sums(pseq_species_noplasmid))

species_p <- sample_sums(pseq_species_plasmid)/sample_sums(pseq_species)
species_p_df <- data.frame(SamplePlasmidProportion=species_p)
ggplot(species_p_df, aes(SamplePlasmidProportion)) + geom_histogram()
ggsave("species_plasmid_histogram.png")

species_df_m <- melt(species_df, value.name="Count")

ggplot(species_df_m, aes(x=Sample, y=Count, fill=variable)) + geom_bar(position="stack", stat="identity") + ggtitle("Total plasmid counts vs not-plasmid per sample")

ggsave("plasmid_species.png", width=15, height=5)

phylum_level_phyloseq_noplasmid_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_phylum_ok_drop50k_2018-12-21_nop.rds"
genus_level_phyloseq_noplasmid_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_genus_all_drop50k_2018-12-21_nop.rds"
species_level_phyloseq_noplasmid_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_species_all_drop50k_2018-12-21_nop.rds"

saveRDS(pseq_species_noplasmid, file=species_level_phyloseq_noplasmid_path)
saveRDS(pseq_genus_noplasmid, file=genus_level_phyloseq_noplasmid_path)
saveRDS(pseq_phylum_noplasmid, file=phylum_level_phyloseq_noplasmid_path)

# add hfc index
hfc <- read.csv("/csc/fr_metagenome/nutrition/hfc_score_barcode.csv")


hfc_adder <- function(pseq, hfc) {
    rownames(hfc) <- hfc[["Barcode"]]

    pheno <- sample_data(pseq_genus_noplasmid)

    hfc_df <- data.frame("HFC_score"=hfc[rownames(pheno), "HFC_score"])
    pheno2 <- cbind(pheno, hfc_df)
    pheno2 <- sample_data(pheno2)

    sample_data(pseq) <- pheno2

    pseq

}

phylum_level_phyloseq_noplasmid_hfc_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_phylum_ok_drop50k_2018-12-21_nop_hfc.rds"
genus_level_phyloseq_noplasmid_hfc_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_genus_all_drop50k_2018-12-21_nop_hfc.rds"
species_level_phyloseq_noplasmid_hfc_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_species_all_drop50k_2018-12-21_nop_hfc.rds"

pseq_genus_noplasmid <- hfc_adder(pseq_genus_noplasmid, hfc)
pseq_species_noplasmid <- hfc_adder(pseq_species_noplasmid, hfc)
pseq_phylum_noplasmid <- hfc_adder(pseq_phylum_noplasmid, hfc)

saveRDS(pseq_species_noplasmid, file=species_level_phyloseq_noplasmid_hfc_path)
saveRDS(pseq_genus_noplasmid, file=genus_level_phyloseq_noplasmid_hfc_path)
saveRDS(pseq_phylum_noplasmid, file=phylum_level_phyloseq_noplasmid_hfc_path)


## Remove antibiotic users and NAs from non-plasmid

# investigate 
sum(meta(pseq_species_noplasmid)[["BL_USE_RX_J01"]], na.rm=TRUE) # 905 cases, 115 NAs

bl_use_rxj01_or_na <- meta(pseq_species_noplasmid)[["BL_USE_RX_J01"]] == 1 | is.na(meta(pseq_species_noplasmid)[["BL_USE_RX_J01"]])

pseq_species_noplasmid_norxj01 <- prune_samples(!bl_use_rxj01_or_na, pseq_species_noplasmid)
pseq_genus_noplasmid_norxj01 <- prune_samples(!bl_use_rxj01_or_na, pseq_genus_noplasmid)
pseq_phylum_noplasmid_norxj01 <- prune_samples(!bl_use_rxj01_or_na, pseq_phylum_noplasmid)

phylum_level_phyloseq_noplasmid_norxj01_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_phylum_ok_drop50k_2018-12-21_nop_norxj01.rds"
genus_level_phyloseq_noplasmid_norxj01_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_genus_all_drop50k_2018-12-21_nop_norxj01.rds"
species_level_phyloseq_noplasmid_norxj01_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_species_all_drop50k_2018-12-21_nop_norxj01.rds"

saveRDS(pseq_species_noplasmid_norxj01, file=species_level_phyloseq_noplasmid_norxj01_path)
saveRDS(pseq_genus_noplasmid_norxj01, file=genus_level_phyloseq_noplasmid_norxj01_path)
saveRDS(pseq_phylum_noplasmid_norxj01, file=phylum_level_phyloseq_noplasmid_norxj01_path)

# add hfc index

phylum_level_phyloseq_noplasmid_norxj01_hfc_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_phylum_ok_drop50k_2018-12-21_nop_norxj01_hfc.rds"
genus_level_phyloseq_noplasmid_norxj01_hfc_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_genus_all_drop50k_2018-12-21_nop_norxj01_hfc.rds"
species_level_phyloseq_noplasmid_norxj01_hfc_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/phfinrisk_species_all_drop50k_2018-12-21_nop_norxj01_hfc.rds"


pseq_genus_noplasmid_norxj01 <- hfc_adder(pseq_genus_noplasmid_norxj01, hfc)
pseq_species_noplasmid_norxj01 <- hfc_adder(pseq_species_noplasmid_norxj01, hfc)
pseq_phylum_noplasmid_norxj01 <- hfc_adder(pseq_phylum_noplasmid_norxj01, hfc)

saveRDS(pseq_species_noplasmid_norxj01, file=species_level_phyloseq_noplasmid_norxj01_hfc_path)
saveRDS(pseq_genus_noplasmid_norxj01, file=genus_level_phyloseq_noplasmid_norxj01_hfc_path)
saveRDS(pseq_phylum_noplasmid_norxj01, file=phylum_level_phyloseq_noplasmid_norxj01_hfc_path)


## also remove same people from ko file.

functinal_activities_orig_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/function_activity.rds"
functinal_activities_norxj01_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/data_aaro/function_activity_norxj01.rds"

functional_items <- readRDS(functinal_activities_orig_path)

modules <- functional_items[["modules"]]
pathways <- functional_items[["pathways"]]
ko <- functional_items[["ko"]]

# p x n matrices

modules_norxj01 <- modules[, !bl_use_rxj01_or_na] 
pathways_norxj01 <- pathways[, !bl_use_rxj01_or_na]
ko_norxj01 <- ko[, !bl_use_rxj01_or_na]

functional_itesm_norxj01 <- list(modules=modules_norxj01,
                                 pathways=pathways_norxj01,
                                 ko=ko_norxj01)

saveRDS(functional_itesm_norxj01, file=functinal_activities_norxj01_path)