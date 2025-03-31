# additional check for Healthy Diet Score

library(phyloseq)
library(microbiome)
library(dplyr)

# check that HFC score reflects well PCA variance of diet variance.
# HFC score is composed after transforming diet questionnaire answer scales to "mount per time unit" scale
# for simplicity, in this check these transforms are not done.

pseq <- readRDS(genus_level_phyloseq_path)

# PCA of species

pseq_clr <- microbiome::transform(microbiome::transform(pseq, "compositional"), "clr")

species_clr_abundances <- pseq_clr %>%
  abundances %>%
  t 

clr_pca <- prcomp(species_clr_abundances)

# plot HFC scores on PCA.

clr_pca_df <- data.frame(pca1=clr_pca$x[,1], pca2=clr_pca$x[,2], pca3=clr_pca$x[,3])
clr_pca_df$HFC_score <- meta(pseq_clr)$HFC_score


basesize <- 10

p1 <- ggplot(data=clr_pca_df,
          aes(x=pca1, y=pca2)) +
        geom_point(aes(color=HFC_score), alpha = 1, size = floor(basesize/8)) +   
        labs(x = "PC1", y = "PC2") +
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

ggsave("~/dietpca1.png", p1)


p2 <- ggplot(data=clr_pca_df,
          aes(x=pca1, y=pca3)) +
        geom_point(aes(color=HFC_score), alpha = 1, size = floor(basesize/8)) +   
        labs(x = "PC1", y = "PC3") +
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

ggsave("~/dietpca2.png", p2)

p3 <- ggplot(data=clr_pca_df,
          aes(x=pca2, y=pca3)) +
        geom_point(aes(color=HFC_score), alpha = 1, size = floor(basesize/8)) +   
        labs(x = "PC2", y = "PC3") +
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

ggsave("~/dietpca3.png", p3)


# PCA of diet variables, not derived
# variable names starting FR02, KY100

fr02_dvars <- sample_variables(pseq)[grep("FR02_1", sample_variables(pseq))]
ky100_dvars <- sample_variables(pseq)[grep("KY100", sample_variables(pseq))]


diet_df <- meta(pseq)[,c(fr02_dvars,ky100_dvars)]
diet_df_num <- sapply(diet_df, function(f) {
    if (is.factor(f))
        as.numeric(levels(f))[f]
    else
        as.numeric(f)
}) %>% as.data.frame
diet_pca <- prcomp(diet_df_num[complete.cases(diet_df_num),])

diet_pca_df <- data.frame(pca1=diet_pca$x[,1], pca2=diet_pca$x[,2],
    HFC_score=meta(pseq)[complete.cases(diet_df_num),"HFC_score"])

pd1 <- ggplot(data=diet_pca_df,
          aes(x=pca1, y=pca2)) +
        geom_point(aes(color=HFC_score), alpha = 1, size = floor(basesize/8)) +   
        labs(x = "diet variables PC1", y = "diet variables PC2") +
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
              ) + ggtitle("HFC score in PCA of FR02 diet variables")
ggsave("~/dietpca_dietvars.png", pd1)
