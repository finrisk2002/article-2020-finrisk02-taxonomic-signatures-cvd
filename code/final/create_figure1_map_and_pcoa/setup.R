analysis_date <- paste0(Sys.Date())

pltdir <- paste0('output/figures/final/figure1_map_stats_pcoa/', analysis_date, '/')
if (!dir.exists(pltdir)) {
  dir.create(pltdir, recursive=TRUE)
}

source("code/final/main_setup.R")
