## Run the entire pipeline

# set_up paths.
source("code/final/main_setup.R")

## PCoA analyses
source("code/final/ordination_analysis/braycurtis_pcoa_all.R")
# see code/final/ordination_analysis/readme.md for how to compute BrayCurtis efficiently on Atlas

## PCA
source("code/final/PCA/pca_all.R")

## Survival analysis
source("code/final/survival_analysis/survival_main.R")

## Network analysis
source("code/final/network_analysis/compute_network_all.R")

## Distribution tables for supplements
source("code/final/distribution_tables/compute_all.R")
