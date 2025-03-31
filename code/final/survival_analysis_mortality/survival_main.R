## Packages ************************************************************* ####

survival_analysis_packages <- c("microbiome", "survival", "tidyverse", "magrittr", "RColorBrewer", 
                                "knitr", "smoothHR", "reshape2", "cowplot", "ggsci", "openxlsx", 
                                "randomForestSRC", "survivalROC", "ggsignif", "cowplot", 
                                "ggfortify", "grid", "gridExtra", "forestplot")

# Install missing packages
missing_packages <- survival_analysis_packages[!(survival_analysis_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages) != 0) install.packages(missing_packages)

# Load packages
for(p in survival_analysis_packages) library(p, character.only = TRUE)

## Settings ************************************************************* ####

# source main_setup for paths
if (!exists("genus_level_phyloseq_path"))
	source("code/final/main_setup.R")

# pltdir setup

source("code/final/survival_analysis/set_up.R")

# set seed
seed <- 11235

# Control baseline age, bmi and sex
covariates <- c("BL_AGE","BMI", "MEN",
                "CURR_SMOKE", "PREVAL_DIAB", "BL_USE_RX_L", 
                "SYSTM", "BP_TREAT")


alpha_level <- 0.05                  # FDR level
normalize <- TRUE                    # Normalize variables before Cox regression
splines <- TRUE                      # Use splines in Cox regression
status <- "DEATH"                    # Cox regression status
time_to_event <- "DEATH_AGEDIFF"     # Cox regression time to event variable

# Load functions
source("code/final/R_functions.R")

## Cox regression  ****************************************************** ####

# Individual genera
source("code/final/survival_analysis/genus_cox.R")

# Subnet abundances
source("code/final/survival_analysis/subnet_cox.R")

if (!exists("run_functional_cox") || run_functional_cox == TRUE) {
    # Functional groups
    source("code/final/survival_analysis/functional_cox.R")
} else {
  message("run_functional_cox == FALSE, skipping functional cox\n")
}

# PCA 1-3, Shannon, Observed
source("code/final/survival_analysis/PCA_cox.R")

## Survival Random Forest *********************************************** ####

# Random survival forest for core genera, covariates and core + covariates
source("code/final/survival_analysis/survival_random_forest.R")

source("code/final/survival_analysis/survival_random_forest_cross_validation.R")


## Figures ************************************************************** ####
if (!exists("create_figures") || create_figures == TRUE) {

# Fig. 2: PCs as predictors of death
source("code/final/survival_analysis/create_figure_pca_mortality_association.R")

# Ext. Fig. 2: PC3 in eastern/western Finland
source("code/final/survival_analysis/create_figure_PC3_east_west_mortality.R")

# Ext. Fig. 3: PC 1-3 driver species
source("code/final/survival_analysis/create_figure_PC_drivers.R")

# Ext. Fig. 4: PC 3 drivers in eastern/western Finland
source("code/final/survival_analysis/create_figure_PC3_drivers_east_west.R")

# Ext. Fig. 5: Individual genera and core importance
source("code/final/survival_analysis/create_figure_genus_panel_srf_importance.R")

}