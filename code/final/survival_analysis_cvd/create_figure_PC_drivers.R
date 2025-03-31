
## Data
species_pseq <- readRDS(file = species_level_phyloseq_path)

pca_clr <- get_pca(species_pseq)

# PC driver taxa
PC_drivers <- lapply(1:3, function(x) {
  
  df <- pca_clr$rotation[, c(x, x+1)] %>%
    as.data.frame() %>% 
    rownames_to_column(var="genus")
  
  genera <- df[, 1:2] %>% set_colnames(c("genus", "PC")) %>% 
    arrange(desc(abs(PC))) %>% 
    mutate(PC = round(PC, 3))
  
}) %>% set_names(paste0("PC", 1:3))

# Combine top20 drivers in a data frame
PC_drivers_df <- lapply(1:length(PC_drivers), function(i){
  n = 20
  
  x <- PC_drivers[[i]] %>%
    arrange(desc(abs(PC))) %>% 
    head(n = n) %>% 
    arrange(desc(PC))
  
  x$genus <- x$genus %>%
    gsub(" \\(Bacteria\\)", "",.) %>% 
    gsub("_", " ", .) %>% 
    gsub("BacteriaPlasmid", "Bacteria Plasmid", .)
    

  
  
  x %>% set_colnames(paste0(c("Genus", "Loading"), " (PC", i,")"))
  
}) %>%
  do.call(cbind, .)


# csv and excel
write_csv(PC_drivers_df, path = "code/survival/shogun_results/PC_drivers.csv")
write.xlsx(PC_drivers_df,
           file = "code/survival/shogun_results/PC_drivers.xlsx", colWidths = "auto")


## Drivers panel ************************
genus_PC_plots <- lapply(1:length(PC_drivers), function(i) {
  
  df <- PC_drivers_df[, (2*i - 1):(2*i)] %>% 
    set_colnames(c("genus", "PC")) %>% 
    arrange(desc(abs(PC)))
  
  
  # Invert axes if more negatives than positives
  # if(mean(df$PC < 0) > .5) {
  #   df$PC <- -1*df$PC
  # }
  
  
  df[, 1] <- as.factor(df[, 1])
  
  
  
  
  
  df %>%
    ggplot(aes(x=factor(genus, levels = rev(df$genus)), y=PC)) +
    geom_col(color = "black") +
    labs(x="", y=paste0("PC", i, " Loading"), title = "") + 
    theme_classic(20) +
    theme(axis.text.x = element_text(angle=45, hjust = 1)) +
    coord_flip()
  
  
})


p <- plot_grid(genus_PC_plots[[1]], 
          genus_PC_plots[[2]], 
          genus_PC_plots[[3]], 
          nrow = 1)


# Save figure
if (!dir.exists(paste0(pltdir, "figure_PC_drivers/"))) {
  dir.create(paste0(pltdir, "figure_PC_drivers/"), recursive=TRUE)
}

png(filename = paste0(pltdir,"figure_PC_drivers/PC_drivers.png"), 
    res = 300, width = 25, height = 10, units = "in")
p
dev.off()