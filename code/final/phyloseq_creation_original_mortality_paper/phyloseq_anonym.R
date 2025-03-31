# this script produces anonymized phyloseqs and writes them to particular location
# * further anonymizes samples by re-id and reorder.
# * adds random noise + categorization to some variables for further anonymization
# * strips our working phyloseqs from phenotype data that the Institute of Health and Welfare does not allow to be shared publicly
# Aaro Salosensaari, May 2020.


library(phyloseq)
library(microbiome)

# set input paths
source("code/final/main_setup.R")

# output file paths, adjust if necessary

phylum_level_phyloseq_public_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/fr02_mortality_anonymized_data/phyloseq/pseq_phylum.rds"
genus_level_phyloseq_public_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/fr02_mortality_anonymized_data/phyloseq/pseq_genus.rds"
species_level_phyloseq_public_path <- "/csc/fr_metagenome/microbiome_scratch/scratch/fr02_mortality_anonymized_data/phyloseq/pseq_species.rds"

outdir_shogun <- "/csc/fr_metagenome/microbiome_scratch/scratch/fr02_mortality_anonymized_data/shogun_tables/"
outdir_pheno <- "/csc/fr_metagenome/microbiome_scratch/scratch/fr02_mortality_anonymized_data/phenotype/"

variables <- c("EAST", "MEN", "BL_AGE","BMI", "DEATH", "DEATH_AGEDIFF")
variables_out <- c("EAST", "MEN", "BL_AGE_cat", "BMI_cat", "DEATH", "DEATH_AGEDIFF_anon")

# shogun phyloseqs, remove other variables
my_helper <- function(pseq_path, variables) {
	pseq <- readRDS(pseq_path)
	variables_in_pseq <- variables[variables %in% sample_variables(pseq)]
	sample_data(pseq) <- sample_data(meta(pseq)[,variables_in_pseq])
	pseq
}


pseq_phylum <- my_helper(phylum_level_phyloseq_path, variables)
pseq_genus  <- my_helper(genus_level_phyloseq_path, variables)
pseq_species <- my_helper(species_level_phyloseq_path, variables)

# extract OTU (phylum, genus, species) and pheno tables
otu_phylum <- abundances(pseq_phylum)
otu_genus <- abundances(pseq_genus)
otu_species <- abundances(pseq_species)

pheno <- meta(pseq_species)


# in new pheno, cut BMI range. also remove NA BMIs (2 items)
pheno[,"BMI_cat"] <- cut(pheno[, "BMI"], breaks=c(0, 25, 30, 60),
    labels=c("Lean", "Overweight", "Obese"), right=FALSE)

# remove all NA
pheno <- pheno[!is.na(pheno$BMI_cat), ]

pheno <- pheno[!is.na(pheno$DEATH), ]


# noise to anonymize

# original summary
#> summary(pheno$BL_AGE)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  24.10   39.07   50.60   49.47   59.80   74.24 


pheno[, "BL_AGE_cat"] <- cut(pheno[, "BL_AGE"], breaks=c(24,40,60, 75), labels=c("24_39", "40_59", "60_75"), right=FALSE, include.lowest=TRUE)


# > summary(pheno$DEATH_AGEDIFF)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.03   14.78   14.84   14.22   14.89   14.96     115

# people who die, anonymize death day

pheno[, "DEATH_AGEDIFF_anon"] <- pheno[, "DEATH_AGEDIFF"] 
pheno[pheno$DEATH==1,"DEATH_AGEDIFF_anon"] <- pheno[pheno$DEATH==1,"DEATH_AGEDIFF_anon"] +  rnorm(sum(pheno$DEATH==1, na.rm=T), mean =0, sd=3)
pheno[, "DEATH_AGEDIFF_anon"] <- round(pheno[, "DEATH_AGEDIFF_anon"])
pheno[pheno[, "DEATH_AGEDIFF_anon"] > 15 , "DEATH_AGEDIFF_anon"] <- 15
pheno[pheno[, "DEATH_AGEDIFF_anon"] < 0  , "DEATH_AGEDIFF_anon"] <- 0

# remove too identifiable subgroups
library(data.table)
dt <- as.data.table(pheno)

dta <- dt[, .N, by=c("MEN", "DEATH", "EAST", "BMI_cat", "BL_AGE_cat")]



# filter rare groups, better than BL_AGE filters (pheno2, pheno3)
groups_remove <- dta[dta$N<30, c("MEN", "DEATH", "EAST", "BMI_cat", "BL_AGE_cat")]

pheno_rmv <- pheno

for (i in 1:nrow(groups_remove)) {
    for (j in 1:nrow(pheno)) {
        if (all(pheno[j, c("MEN", "DEATH", "EAST", "BMI_cat", "BL_AGE_cat")] == as.data.frame(groups_remove[i,])))
            pheno_rmv[j,] <- NA
    }
    
}

# rmeove uncomplete cases

pheno <- pheno_rmv[complete.cases(pheno_rmv),]

pheno <- pheno[, variables_out]

otu_phylum <- otu_phylum[, rownames(pheno)]
otu_genus <- otu_genus[, rownames(pheno)]
otu_species <- otu_species[, rownames(pheno)]


## reorder

# find a reordering for sample names
sample_ids <- rownames(pheno)
reord <- sample.int(length(sample_ids))
reord_sids <- sample_ids[reord]

# apply the reordering
otu_phylum <- otu_phylum[, reord_sids]
otu_genus <- otu_genus[, reord_sids]
otu_species <- otu_species[, reord_sids]

pheno <- pheno[reord_sids,]

# set new IDs
oldids <- colnames(otu_genus)
newids <- paste0("s", 1:length(sample_ids))
colnames(otu_phylum) <- newids
colnames(otu_genus) <- newids
colnames(otu_species) <- newids

rownames(pheno) <- newids




# create new phyloseqs
sam <- sample_data(pheno)
sample_names(sam) <- newids

pseq_phylum <- phyloseq(otu_table(otu_phylum, taxa_are_rows=TRUE),
						sam, tax_table(pseq_phylum))

pseq_genus <- phyloseq(otu_table(otu_genus, taxa_are_rows=TRUE),
						sam, tax_table(pseq_genus))

pseq_species <- phyloseq(otu_table(otu_species, taxa_are_rows=TRUE),
						sam, tax_table(pseq_species))

# see main_setup.R
saveRDS(pseq_phylum, phylum_level_phyloseq_public_path)
saveRDS(pseq_genus, genus_level_phyloseq_public_path)
saveRDS(pseq_species, species_level_phyloseq_public_path)


write.csv(pheno, file=paste0(outdir_pheno, "/phenotype.csv"), row.names=TRUE)
write.csv(abundances(pseq_genus), file=paste0(outdir_shogun, "/table_genus.csv"), row.names=TRUE)
write.csv(abundances(pseq_species), file=paste0(outdir_shogun, "/table_species.csv"), row.names=TRUE)
write.csv(abundances(pseq_phylum), file=paste0(outdir_shogun, "/table_phylum.csv"), row.names=TRUE)

file.copy("code/final/phyloseq_creation/data_anonymization_readme.md",
          "/csc/fr_metagenome/microbiome_scratch/scratch/fr02_mortality_anonymized_data/data_anonymization_readme.md")