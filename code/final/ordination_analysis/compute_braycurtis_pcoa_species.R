library(biomformat)
library(phyloseq)
library(ggplot2)
library(microbiome)

# equivalent to vegan + ape procedure.

phfs.pseq <- readRDS(species_level_phyloseq_path)

phfs.pseq.rel <- microbiome::transform(phfs.pseq, 'compositional')

phfs.pseq.rel.pcoa <- ordinate(phfs.pseq.rel, method="PCoA", distance="bray")

saveRDS(phfs.pseq.rel.pcoa, file=paste0(datadir, "phf_species_d50k_bray_pseq_pcoa",
                                 analysis_date, ".rds"))