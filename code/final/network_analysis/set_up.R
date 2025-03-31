library(biomformat)
library(microbiome)
library(phyloseq)
library(data.table)
library(RColorBrewer)
library(dplyr)
packageVersion('phyloseq')
# library(finriskmetagcommon) # DEPRECATE
library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(grid)
library(gtable)

library(survival)
library(smoothHR)
library(tibble)
library(magrittr)

library(scales)

library(network)
library(sna)
library(GGally)
library(SpiecEasi)


analysis_date <- paste0(Sys.Date())


if (!exists("analysis_ident"))
    source("code/final/main_setup.R")

pltdir <- paste0('output/figures/final/figure_subnets_shogun/', analysis_date, '/', analysis_ident, '/')
if (!dir.exists(pltdir)) {
  dir.create(pltdir, recursive=TRUE)
}


# for tables
outdir <- paste0('output/tables/final/spieceasi_tables/', analysis_date, '/', analysis_ident, '/')
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive=TRUE)
}
