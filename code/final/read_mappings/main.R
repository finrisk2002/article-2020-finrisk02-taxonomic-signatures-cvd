library(rmarkdown)
library(readxl)

# Now use only these samples (the ones from our manuscript)
# Only keep the samples that are in FR02
#include.samples0 <- as.character(unname(unlist(read.table("FINRISK_sample_ids"))))
include.samples <- as.character(unname(unlist(read_excel("shogun_sample_ids.xlsx"))))

source("read_bam.R"); save.image()

source("analyse.R"); save.image()

source("stats.R"); save.image()

source("report.R");

