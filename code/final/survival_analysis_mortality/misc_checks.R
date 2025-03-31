
pseq <- readRDS(file = genus_level_phyloseq_path)
species_pseq <- readRDS(file = species_level_phyloseq_path)
phylum_pseq <- readRDS(file = phylum_level_phyloseq_path)

full_pseq <- readRDS(file = genus_level_NO_drop50k_phyloseq_path)



# pseq <- readRDS(file = genus_level_phyloseq_path)
# species_pseq <- readRDS(file = species_level_phyloseq_path)

## Number of taxa with different thresholds *********** ####

# Genus 20%, 50% and 90%
pseq_compositional <- microbiome::transform(pseq, "compositional")
phylum_compositional <- microbiome::transform(phylum_pseq, "compositional")
species_compositional <- microbiome::transform(species_pseq, "compositional")

thresholds <- c(20, 50, 90)

n_species <- lapply(thresholds, function(prop) {
  
  core(species_compositional,
       detection = 0,
       prevalence = prop/100) %>% 
    ntaxa
  
}) %>% 
  unlist() %>% set_names(paste0(thresholds, "%"))

n_genus <- lapply(thresholds, function(prop) {
  
  core(pseq_compositional,
       detection = 0,
       prevalence = prop/100) %>% 
    ntaxa
  
}) %>% 
  unlist() %>% set_names(paste0(thresholds, "%"))

n_phylum <- lapply(thresholds, function(prop) {
  
  core(phylum_compositional,
       detection = 0,
       prevalence = prop/100) %>% 
    ntaxa
  
}) %>%
  unlist() %>%
  set_names(paste0(thresholds, "%"))


thresholds2 <- c(0, 1)

n_species2 <- lapply(thresholds2, function(prop) {
  
  core(species_compositional,
       detection = 0.1/100,
       prevalence = prop/100) %>% 
    ntaxa
  
}) %>% 
  unlist() %>% set_names(paste0(thresholds2, "%"))

n_genus2 <- lapply(thresholds2, function(prop) {
  
  core(pseq_compositional,
       detection = 0.1/100,
       prevalence = prop/100) %>% 
    ntaxa
  
}) %>% 
  unlist() %>% set_names(paste0(thresholds2, "%"))

n_phylum2 <- lapply(thresholds2, function(prop) {
  
  core(phylum_compositional,
       detection = 0.1/100,
       prevalence = prop/100) %>% 
    ntaxa
  
}) %>%
  unlist() %>%
  set_names(paste0(thresholds2, "%"))



## EAST/WEST common genera **************************** ####

n_common_genera <- intersect(taxa_names(core(subset_samples(pseq, EAST == "EAST"),
                                             detection = 0/100, 
                                             prevalence = 0/100)), 
                             taxa_names(core(subset_samples(pseq, EAST == "WEST"),
                                             detection = 0/100, 
                                             prevalence = 0/100))) %>%
  length

## 95% core ******************************************* ####
core95 <- core(
  pseq_compositional,
  detection = 0, prevalence = 95/100, 
)

compositional_df <- abundances(pseq_compositional)

core95_abundances <- lapply(1:ncol(compositional_df), function(i) {
  
  compositional_df[taxa_names(core95), i] %>% 
    sum
  
}) %>% unlist

median(core95_abundances)

## Read histogram ************************************* ####

counts <- colSums(abundances(full_pseq), na.rm = T) %>% 
  data.frame(counts = .)


counts %>% 
  ggplot() + 
  geom_histogram(aes(x = log10(counts)), 
                 fill = "white",
                 color = "black",
                 bins = 200) +
  geom_vline(xintercept = log10(5e4), linetype = "dashed", color = "red") +
  labs(x = expression(Log[10]~("Read Count")), y = "Frequency")
# coord_cartesian(xlim = c(0, 1e7)) +
# scale_x_continuous(breaks = c(1e6, 5e6, 1e7))

## Infection death ************************************ ####

icd <- c(paste0("A0", 47:99),
         # paste0("A", 47:99),
         paste0("A", 100:499), 
         "A159", 
         "B250", "B377", "B59", 
         "G051", "G008", 
         "I330", "I520", 
         paste0("J", 13:16), paste0("J", 18:20), 
         "J40", "J69", "J171", "J172", 
         paste0("K", 65:81), 
         "L031",
         "M861", 
         "N10", "N390", 
         "T874")



meta(pseq)[ , "K_TPKS"] %>%
  table %>%
  as.data.frame() %>% 
  set_colnames(c("code", "Freq")) %>% 
  filter(code %in% icd)
## Centrifuge vs. Shogun PC4/3 Spearman *************** ####

centrifuge_pca <- readRDS(file = centrifuge_genus_level_phyloseq_path)
shogun_pca <- readRDS(file = genus_level_phyloseq_path)


# Check that PC4 in Enterobacteriaceae driven
centrifuge_drivers <- centrifuge_pca$rotation %>%
  as.data.frame() %>% 
  rownames_to_column("driver") %>%
  select(1:5) %>% 
  arrange(desc(abs(PC4))) 


centrifuge_PC4 <- centrifuge_pca$x %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  select(ID, centrifuge_PC4 = PC4)

shogun_PC3 <- shogun_pca$x %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  select(ID, shogun_PC3 = PC3)


# Combine
death_PC_df <- full_join(centrifuge_PC4, shogun_PC3, by = "ID")

cor(-death_PC_df$centrifuge_PC4,
    death_PC_df$shogun_PC3,
    method = "spearman",
    use = "complete.obs")



