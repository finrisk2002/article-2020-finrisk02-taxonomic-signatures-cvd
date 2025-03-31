## ----setup, include=FALSE------------------------------------------------

library(rmarkdown)
library(biomformat)
library(phyloseq)
library(ggplot2)
library(microbiome)
library(RColorBrewer)
library(grid)
#library(finriskmetagcommon)
library(dplyr)

## ------------------------------------------------------------------------

source('pcoa_functions.R')

phf <- phfinrisk_species
phf.rel <- microbiome::transform(phf, 'compositional')

## ------------------------------------------------------------------------

tt.m <- as(tax_table(phf.rel), "matrix")
phyla.names <- na.omit(unique(tt.m[, "Phylum"]))

# how much the abundance of each phyla correlates with the first two PCoA axis?
# also, we skip Plasmid phyla
phyla.prop.abd <- lapply(phyla.names,
                         function (x) {
           phyl.taxa <- which(tt.m[, "Phylum"] == x & !grepl("Plasmid", tt.m[, "Kingdom"])) %>%
	     as.integer
           # abundances of taxa
           phyl.abud <- apply(otu_table(phf.rel)[phyl.taxa, ], 2, function (x) sum(x))
                         })
		 
phyla.df <- do.call(cbind.data.frame, phyla.prop.abd)
ab.names <- lapply(phyla.names, function (x) paste0(x,'.Abundance') )
colnames(phyla.df) <- ab.names
phyla.df[, 1:length(phyla.names)] <- sapply(phyla.df[, 1:length(phyla.names)], as.numeric)
phyla.df$pcoa1 <- phfs.rel.ord.pcoa$vectors[,1]
phyla.df$pcoa2 <- phfs.rel.ord.pcoa$vectors[,2]

# Genus-level plots
prevo.taxa <- tt.m[,"Genus"] == "Prevotella"
prevotella.abud <- apply(otu_table(phf.rel)[prevo.taxa, ], 2, function (x) sum(x, na.rm=T))
phyla.df$Prevotella.Abundance <- prevotella.abud %>% unlist



# Calculate species loadings on PCoA.
skip <- TRUE
if (!skip) {

  o <- abundances(transform(phfinrisk_genus, "compositional"))
  # Takes very long time to calculate
  # p <- plot_landscape(transform(phfinrisk_genus, "compositional"), distance = "bray", method = "PCoA")
  # ggsave('phfs.rel.ord.pcoa.prevotella.png', width=5, height=5, dpi=700)

  # Sample similarities
  d <- vegdist(t(o))
  
  # PCoA and percent of explained variation w.r.t. the full dimensional data
  #pc.bray2 <- cmdscale(d,k=2)
  # BrayCurtis is not a metric. Hence use add argument to avoid negative eigenvalues.
  pc.bray <- cmdscale(d, k=(ncol(o)-1), eig = TRUE, add = TRUE)
  eig2 <- eigenvals(pc.bray)
  perc <- eig2 / sum(eig2)

  # Significant PCoA dimensions
  #library(PCPS)
  #library(parallel)  
  #signif <- pcoa.sig(t(o), method = "bray", axis = 10, iterations = 100, parallel = detectCores())

  # Standard method: each species as a weighted average of the positions
  #  of the samples in which it is present. The weights are the relative
  #  abundances of that species in the samples.
  #o.norm <- sweep(t(o), 2, colSums(o), '/')
  #wa <- t(o.norm) %*% pc.bray
  #par(mar = c(2, 15, 1,1 )); s <- sort(wa[,1]); barplot(c(head(s), tail(s)), horiz = T, las = 1)
  #par(mar = c(2, 15, 1,1 )); s <- sort(wa[,2]); barplot(c(head(s), tail(s)), horiz = T, las = 1)

}

#- -------------------------------------

# PCoA in 3D
skip <- TRUE
if (!skip) {
  pc3 <- cmdscale(d.bray,k=3)
  # Calculate species loadings on PCoA.
  # Standard method: each species as a weighted average of the positions
  #  of the samples in which it is present. The weights are the relative
  #  abundances of that species in the samples.
  w3 <- t(o.norm) %*% pc3
  par(mar = c(2, 15, 1,1 )); s <- sort(w3[,1]); barplot(c(head(s), tail(s)), horiz = T, las = 1)
  par(mar = c(2, 15, 1,1 )); s <- sort(w3[,2]); barplot(c(head(s), tail(s)), horiz = T, las = 1)
  par(mar = c(2, 15, 1,1 )); s <- sort(w3[,3]); barplot(c(head(s), tail(s)), horiz = T, las = 1)
}
