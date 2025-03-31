
## Load data

# Load phyloseq
pseq <- readRDS(file = genus_level_phyloseq_path)
species_pseq <- readRDS(file = species_level_phyloseq_path)

# Add PCA components and diversity measures
pseq <- add_pca_components_to_meta_data(pseq, species_pseq, n_axes = 3)


## Death ~ PC1-3 panels


PC1_hazard <- plot_HR(fit = coxph(Surv(meta(pseq)$DEATH_AGEDIFF, meta(pseq)$DEATH) ~
                                    BL_AGE + BMI + MEN + CURR_SMOKE + PREVAL_DIAB + BL_USE_RX_L +
                                    SYSTM + BP_TREAT + PC1,
                                  data=meta(pseq),
                                  x=TRUE),
                      pred = "PC1",
                      data=meta(pseq),
                      title = FALSE) +
  labs(x="PC1") +
  theme_classic(20) +
  coord_cartesian(ylim=c(-1,2))

PC2_hazard <- plot_HR(fit = coxph(Surv(meta(pseq)$DEATH_AGEDIFF, meta(pseq)$DEATH) ~
                                    BL_AGE + BMI + MEN + CURR_SMOKE + PREVAL_DIAB + BL_USE_RX_L +
                                    SYSTM + BP_TREAT + PC2,
                                  data=meta(pseq),
                                  x=TRUE),
                      pred = "PC2",
                      data=meta(pseq),
                      title = FALSE) +
  labs(x="PC2", y= "") +
  theme_classic(20) +
  coord_cartesian(ylim=c(-1,2)) + 
  scale_y_continuous(labels = NULL) 

PC3_hazard <- plot_HR(fit = coxph(Surv(meta(pseq)$DEATH_AGEDIFF, meta(pseq)$DEATH) ~
                                    BL_AGE + BMI + MEN + CURR_SMOKE + PREVAL_DIAB + BL_USE_RX_L +
                                    SYSTM + BP_TREAT + PC3,
                                  data=meta(pseq),
                                  x=TRUE),
                      pred = "PC3",
                      data=meta(pseq),
                      title = FALSE) +
  labs(x="PC3", y= "") +
  theme_classic(20) +
  coord_cartesian(ylim=c(-1,2)) + 
  scale_y_continuous(labels = NULL) 



p <- plot_grid(PC1_hazard, 
               PC2_hazard, 
               PC3_hazard, 
               nrow = 1, rel_widths = c(1.1, 1, 1))



# Save
if (!dir.exists(paste0(pltdir, "figure_mortality_vs_PCs/"))) {
  dir.create(paste0(pltdir, "figure_mortality_vs_PCs/"), recursive=TRUE)
}
png(filename = paste0(pltdir, "figure_mortality_vs_PCs/mortality_vs_PCs.png"), 
    res = 300, width = 10, height = 3.5, units = "in")
p
dev.off()


