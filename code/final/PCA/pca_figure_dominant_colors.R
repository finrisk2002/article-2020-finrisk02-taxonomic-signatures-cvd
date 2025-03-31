# Aaro 2020-01-22
pseq.genus <- readRDS(genus_level_phyloseq_path)

species_pseq <- readRDS(species_level_phyloseq_path)


# Compositional
species_pseq_compositional <- microbiome::transform(species_pseq, "compositional")

species_pseq_clr <- microbiome::transform(species_pseq_compositional, "clr")

# CLR abundances
species_clr_abundances <- species_pseq_clr %>%
  abundances %>%
  t 


# PCA
species_clr_pca <- prcomp(species_clr_abundances)

### ...

pseq.rel <- species_pseq_compositional

pseq.genus.rel <- microbiome::transform(pseq.genus, "compositional")

# find and plot the dominant phyla in the species-level 2-d PCoA

tt.m <- as(tax_table(pseq.rel), "matrix")

phyla.names <- na.omit(unique(tt.m[, "Phylum"]))

# obtain phylum-level proportional abundances
# also, we skip Plasmid phyla
phyla.prop.abd <- lapply(phyla.names,
                         function (x) {
                           phyl.taxa <- which(tt.m[, "Phylum"] == x &
                            !grepl("Plasmid", tt.m[, "Kingdom"])) %>% as.integer
                           # abundances of taxa
                           phyl.abud <- apply(otu_table(pseq.rel)[phyl.taxa, ], 2, function (x) sum(x))
                         })

phyla.df <- do.call(cbind.data.frame, phyla.prop.abd)
ab.names <- lapply(phyla.names, function (x) paste0(x,'.Abundance') )
colnames(phyla.df) <- ab.names
phyla.df[, 1:length(phyla.names)] <- sapply(phyla.df[, 1:length(phyla.names)], as.numeric)

df <- phyla.df
df$dominant.genera <- taxa(pseq.genus)[apply(abundances(pseq.genus.rel), 2, function (x) {which.max(x)})]
df$dominant.genera <- gsub("\\(Bacteria\\)", "", df$dominant.genera)
df$dominant.genera[!df$dominant.genera %in% names(rev(sort(table(df$dominant.genera))))[1:6]] <- "Other"
df$dominant.genera <- factor(df$dominant.genera, levels = names(sort(table(df$dominant.genera))))
df <- bind_cols(df, meta(pseq.rel))

#pca <- princomp(t(abundances(transform(phfinrisk_genus, "clr"))))
#df$pca1 <- pca$scores[, 1]
#df$pca2 <- pca$scores[, 2]

pca <- species_clr_pca
df$pca1 <- pca$x[,1]
df$pca2 <- pca$x[,2]
df$pca3 <- pca$x[,3]

# skip plot creation if create_figures == FALSE. default behavior: create figures

if (!exists("create_figures") || create_figures == TRUE) {

basesize <- 100

p <- ggplot(data=subset(df, dominant.genera %in% rev(levels(df$dominant.genera))[1:length(unique(df$dominant.genera))]) %>%
        arrange(desc(dominant.genera)),
          aes(x=pca1, y=pca2)) +
        geom_point(aes(color=dominant.genera), alpha = 1, size = floor(basesize/8)) +   
        labs(x = "PC1", y = "PC2") +
        scale_colour_manual(values = c("black", "blue", "lightblue", "darkgray", "magenta", "darkgreen", "red")) +
        guides(color = guide_legend(title = "Dominant genus", reverse = TRUE, keyheight=0.35, default.unit="inch")) +
        theme_classic(base_family = "", base_size = basesize) +  	
        theme(# legend.position = "left",
              legend.position = c(0.165, 0.22),
              legend.text = element_text(size = 1 * basesize),
              legend.title = element_text(size = 1 * basesize),
              legend.background = element_rect(fill="transparent"),
              axis.text.x = element_text(size = basesize),
              axis.text.y = element_text(size = basesize),
              axis.title.x = element_text(size = basesize),
              axis.title.y = element_text(size = basesize)
              )

figure.pca <- p

print(p)

theme_set(theme_bw(basesize))
ggsave(paste0(pltdir, "pca_dominant_colors.png"), width=40, height=40)
ggsave(paste0(pltdir, "pca_dominant_colors.pdf"), width=40, height=40)


p <- ggplot(data=subset(df, dominant.genera %in% rev(levels(df$dominant.genera))[1:length(unique(df$dominant.genera))]) %>%
        arrange(desc(dominant.genera)),
          aes(x=pca1, y=pca3)) +
        geom_point(aes(color=dominant.genera), alpha = 1, size = floor(basesize/8)) +   
        labs(x = "PC1", y = "PC3") +
        scale_colour_manual(values = c("black", "blue", "lightblue", "darkgray", "magenta", "darkgreen", "red")) +
        guides(color = guide_legend(title = "Dominant genus", reverse = TRUE, keyheight=0.35, default.unit="inch")) +
        theme_classic(base_family = "", base_size = basesize) +  	
        theme(# legend.position = "left",
              legend.position = c(0.165, 0.22),
              legend.text = element_text(size = 1 * basesize),
              legend.title = element_text(size = 1 * basesize),
              legend.background = element_rect(fill="transparent"),
              axis.text.x = element_text(size = basesize),
              axis.text.y = element_text(size = basesize),
              axis.title.x = element_text(size = basesize),
              axis.title.y = element_text(size = basesize)
              )

p

ggsave(paste0(pltdir, "pca_dominant_colors_pc3.png"), width=40, height=40)


}