if (!exists("analysis_ident"))
    source("code/final/main_setup.R")

analysis_date <- paste0(Sys.Date())

pltdir <- paste0('output/figures/final/figure_survival/', analysis_date, '/', analysis_ident, '/')
if (!dir.exists(pltdir)) {
  dir.create(pltdir, recursive=TRUE)
}

out_tabledir <- paste0('output/tables/', analysis_date, '/', analysis_ident, '/')
if (!dir.exists(out_tabledir)) {
  dir.create(out_tabledir, recursive=TRUE)
}