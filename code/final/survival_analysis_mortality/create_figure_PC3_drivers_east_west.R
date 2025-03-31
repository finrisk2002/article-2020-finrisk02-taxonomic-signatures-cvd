
## Load data

# Load phyloseq
pseq <- readRDS(file = genus_level_phyloseq_path)
species_pseq <- readRDS(file = species_level_phyloseq_path)

# Add PCA components and diversity measures
pseq <- add_pca_components_to_meta_data(pseq, species_pseq, n_axes = 3)


## Subset pseqs

# Genus level
east_pseq <- pseq %>% 
  subset_samples(EAST == "EAST")
west_pseq <- pseq %>% 
  subset_samples(EAST == "WEST")

# Species level
east_species_pseq <- species_pseq %>%
  subset_samples(EAST == "EAST")
west_species_pseq <- species_pseq %>%
  subset_samples(EAST == "WEST")



## PCA: separately for both areas
pca_clr_east <- get_pca(east_species_pseq, check_existence = FALSE)
pca_clr_west <- get_pca(west_species_pseq, check_existence = FALSE)


# saveRDS(pca_clr_east, file = "output/raw_output/species_pca_clr_east.Rdata")
# saveRDS(pca_clr_west, file = "output/raw_output/species_pca_clr_west.Rdata")


# EAST
east_PC_drivers <- lapply(1:5, function(x) {
  
  df <- pca_clr_east$rotation[, c(x, x+1)] %>%
    as.data.frame() %>% 
    rownames_to_column(var="genus")
  
  genera <- df[, 1:2] %>% set_colnames(c("genus", "PC")) %>% 
    arrange(desc(abs(PC))) %>% 
    mutate(PC = round(PC, 3))
  
}) %>% set_names(paste0("PC", 1:5))
east_PC_drivers_df <- east_PC_drivers[["PC1"]]
for(i in 2:5) {
  
  east_PC_drivers_df <- full_join(east_PC_drivers_df, east_PC_drivers[[paste0("PC", i)]], "genus")
  
}

east_PC_drivers_df <- east_PC_drivers_df %>% 
  set_colnames(c("genus", paste0("PC", 1:5))) %>% 
  cbind(area = "EAST")


# WEST
west_PC_drivers <- lapply(1:5, function(x) {
  
  df <- pca_clr_west$rotation[, c(x, x+1)] %>%
    as.data.frame() %>% 
    rownames_to_column(var="genus")
  
  genera <- df[, 1:2] %>% set_colnames(c("genus", "PC")) %>% 
    arrange(desc(abs(PC))) %>% 
    mutate(PC = round(PC, 3))
  
}) %>% set_names(paste0("PC", 1:5))
west_PC_drivers_df <- west_PC_drivers[["PC1"]]
for(i in 2:5) {
  
  west_PC_drivers_df <- full_join(west_PC_drivers_df, west_PC_drivers[[paste0("PC", i)]], "genus")
  
}

west_PC_drivers_df <- west_PC_drivers_df %>% 
  set_colnames(c("genus", paste0("PC", 1:5))) %>% 
  cbind(area = "WEST")


## Intersection of top 20 in east and west
common_top_20 <- intersect(west_PC_drivers_df %>% 
                             arrange(desc(abs(PC3))) %>% 
                             head(n = 20) %>%
                             pull(genus), 
                           east_PC_drivers_df %>% 
                             arrange(desc(abs(PC3))) %>% 
                             head(n = 20) %>%
                             pull(genus))

## Plot PC3 bar plots

west_pc3 <- west_PC_drivers_df %>%
  arrange(desc(abs(PC3))) %>% 
  filter(genus %in% common_top_20) %>% 
  select(genus, PC3) %>% 
  mutate(PC3 = abs(PC3),
         Area = "West")

east_pc3 <- east_PC_drivers_df %>%
  arrange(desc(abs(PC3))) %>% 
  filter(genus %in% common_top_20) %>% 
  select(genus, PC3) %>% 
  mutate(PC3 = abs(PC3),
         Area = "East")


east_west_pc3 <- rbind(east_pc3, 
                       west_pc3) %>% 
  mutate(genus = gsub("\\(Bacteria\\)", "", genus)) %>% 
  mutate(genus = gsub("BacteriaPlasmid", "Bacteria Plasmid", genus)) %>% 
  mutate(genus = gsub("_", " ", genus)) %>% 
  mutate(genus = genus %>% as.factor()) %>% 
  mutate(Area = Area %>% factor(levels = c("West", "East")))


p <- east_west_pc3 %>%
  ggplot(aes(x = factor(genus, levels = rev((east_west_pc3 %>% filter(Area == "East"))$genus)), y = PC3, fill = Area)) + 
  geom_col(position = "dodge") +
    coord_flip() +
    scale_fill_manual(values = c("red", "blue")) + 
  labs(x = "Species", y = "PC3 Loading")




# Save

if (!dir.exists(paste0(pltdir, "figure_east_west_PC3_drivers/"))) {
  dir.create(paste0(pltdir, "figure_east_west_PC3_drivers/"), recursive=TRUE)
}
png(filename = paste0(pltdir, "figure_east_west_PC3_drivers/east_west_PC3_drivers.png"), 
    res = 300, width = 8, height = 7, units = "in")
p
dev.off()
  