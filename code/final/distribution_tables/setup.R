analysis_date <- paste0(Sys.Date())

outdir <- paste0('output/tables/final/distribution_tables/', analysis_date, '/', analysis_ident, '/')
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive=TRUE)
}

source("code/final/main_setup.R")