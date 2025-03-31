## Species level PCA

# Set seed so the axes orientations will remain under control
set.seed(12345)

# Load species level phyloseq
species_pseq <- readRDS(species_level_phyloseq_path)

# Compositional transformation
species_pseq_compositional <- microbiome::transform(species_pseq, "compositional")

# CLR transformation
species_pseq_clr <- microbiome::transform(species_pseq_compositional, "clr")

# CLR abundances
species_clr_abundances <- species_pseq_clr %>%
  abundances %>%
  t 


# PCA
clr_pca <- prcomp(species_clr_abundances)

# Flip PC1 and PC3 to get a neat figures 
clr_pca$rotation[, c("PC1", "PC3")] <- -1*clr_pca$rotation[, c("PC1", "PC3")]
clr_pca$x[, c("PC1", "PC3")] <- -1*clr_pca$x[, c("PC1", "PC3")]


# save
saveRDS(clr_pca, file = paste0(datadir,"/shogun_species_pca.rds"))
# load(file = "code/survival/data/shogun_species_pca")
# saveRDS(clr_pca, file = "output/raw_output/species_pca.rds")
