# example: phyloseq creation from csvs


library(phyloseq)
library(microbiome)


outdir_shogun <- "output/fr02_mortality_datapublication/shogun_tables"
pheno_file <- "output/fr02_mortality_datapublication/phenotype/phenotype.csv"

# genus and species phyloseqs

abundances_genus <- read.csv(paste0(outdir_shogun, "/table_genus.csv"),
    header=TRUE, stringsAsFactors=TRUE, check.names=FALSE, row.names=1)

abundances_species <- read.csv(paste0(outdir_shogun, "/table_species.csv"),
    header=TRUE, stringsAsFactors=TRUE, check.names=FALSE, row.names=1)

pheno <- read.csv(pheno_file, header=TRUE, stringsAsFactors=TRUE, check.names=FALSE, row.names=1)


pseq_genus <- phyloseq(otu_table(abundances_genus, taxa_are_rows=TRUE), sample_data(pheno))
pseq_species <- phyloseq(otu_table(abundances_species, taxa_are_rows=TRUE), sample_data(pheno))

pseq_outdir <- "output/fr02_mortality_datapublication/phyloseq"

if (!dir.exists(pseq_outdir)) {
  dir.create(pseq_outdir, recursive=TRUE)
}

saveRDS(pseq_genus, file="output/fr02_mortality_datapublication/phyloseq/pseq_genus.rds")
saveRDS(pseq_species, file="output/fr02_mortality_datapublication/phyloseq/pseq_species.rds")