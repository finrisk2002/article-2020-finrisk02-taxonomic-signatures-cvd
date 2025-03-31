df <- phyla.df
df$dominant.genera <- taxa(phfinrisk_genus)[apply(abundances(transform(phfinrisk_genus, "compositional")), 2, function (x) {which.max(x)})]
df$dominant.genera <- gsub("\\(Bacteria\\)", "", df$dominant.genera)
df$dominant.genera[!df$dominant.genera %in% names(rev(sort(table(df$dominant.genera))))[1:6]] <- "Other"
df$dominant.genera <- factor(df$dominant.genera, levels = names(sort(table(df$dominant.genera))))
df <- bind_cols(df, meta(phf.rel))

# ordination object
ord.pcoa <- phfs.rel.ord.pcoa

p <- ggplot(data=subset(df, dominant.genera %in% rev(levels(df$dominant.genera))[1:length(unique(df$dominant.genera))]) %>%
        arrange(desc(dominant.genera)),
          aes(x=pcoa1, y=pcoa2)) +
        geom_point(aes(color=dominant.genera), alpha = 1, size = 1.8) +   
        labs(x = "PCoA 1", y = "PCoA 2") +
        scale_colour_manual(values = c("black", "blue", "lightblue", "darkgray", "magenta", "darkgreen", "red")) +
        guides(color = guide_legend(title = "", reverse = TRUE, keyheight=0.35, default.unit="inch")) +
        theme_classic(base_family = "", base_size = 34) +  	
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

# extract relative eigenvalues (proportion of variance explained by each axis)
p <- p + xlab(paste("PCoA 1 (", round(100*ord.pcoa$values$Relative_eig[1],1), "%", ")", sep = ""))
p <- p + ylab(paste("PCoA 2 (", round(100*ord.pcoa$values$Relative_eig[2],1), "%", ")", sep = ""))

figure.pcoa <- p

theme_set(theme_bw(10))
png("pcoagenus.png", width = 700, height = 700)
print(figure.pcoa)
dev.off()

