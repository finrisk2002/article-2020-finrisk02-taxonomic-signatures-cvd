analysis_date <- paste0(Sys.Date())

pltdir <- paste0('output/figures/final/figure_ordination/', analysis_date, '/')
if (!dir.exists(pltdir)) {
  dir.create(pltdir, recursive=TRUE)
}

#species_level_phyloseq_path <- "input/data_work/phfinrisk_species_all_drop50k_2018-12-21.RDs"
#datadir <- "output/raw_output/"

source("code/final/main_setup.R")