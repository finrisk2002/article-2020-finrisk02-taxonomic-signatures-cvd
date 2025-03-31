library(biomformat)
library(phyloseq)
library(ggplot2)
library(microbiome)
library(ape)

# Copypaste'd set_up.R here for now. not pretty, but works (easier to run with grun.py)

analysis_date <- paste0(Sys.Date())

pltdir <- paste0('output/figures/final/figure_ordination/', analysis_date, '/')
if (!dir.exists(pltdir)) {
  dir.create(pltdir, recursive=TRUE)
}

#datadir <- "output/raw_output/"

# Copypaste ends.


bray.dist.m <- readRDS(paste0(datadir, "phf_species_d50k_bray_dist_m",
                                 "2020-01-31", ".rds"))

phfs.pseq.bd.pcoa <- pcoa(bray.dist.m)

saveRDS(phfs.pseq.bd.pcoa, file=paste0(datadir, "phf_species_d50k_bray_ape_pcoa",
                                 analysis_date, ".rds"))