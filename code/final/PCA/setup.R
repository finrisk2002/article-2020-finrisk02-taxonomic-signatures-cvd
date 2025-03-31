library(dplyr)
library(phyloseq)
library(microbiome)

analysis_date <- paste0(Sys.Date())

pltdir <- paste0('output/figures/final/figure_pca/', analysis_date, '/', analysis_ident, '/')
if (!dir.exists(pltdir)) {
  dir.create(pltdir, recursive=TRUE)
}

source("code/final/main_setup.R")