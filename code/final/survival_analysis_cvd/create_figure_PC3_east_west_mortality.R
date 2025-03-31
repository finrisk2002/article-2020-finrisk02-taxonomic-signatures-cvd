
## Load data

# Load phyloseq
pseq <- readRDS(file = genus_level_phyloseq_path)
species_pseq <- readRDS(file = species_level_phyloseq_path)

# Add PCA components and diversity measures
pseq <- add_pca_components_to_meta_data(pseq, species_pseq, n_axes = 3)


# Subset pseqs
east_pseq <- pseq %>% 
  subset_samples(EAST == "EAST")
west_pseq <- pseq %>% 
  subset_samples(EAST == "WEST")


east_species_pseq <- species_pseq %>%
  subset_samples(EAST == "EAST")

west_species_pseq <- species_pseq %>%
  subset_samples(EAST == "WEST")

## PCA: separately for both areas

east_pseq <- add_pca_components_to_meta_data(east_pseq, east_species_pseq, n_axes = 3, check_existence = FALSE)
west_pseq <- add_pca_components_to_meta_data(west_pseq, west_species_pseq, n_axes = 3, check_existence = FALSE)



## Plot

east_p <- plot_HR(fit=coxph(Surv(meta(east_pseq)$DEATH_AGEDIFF, meta(east_pseq)$DEATH) ~
                              BL_AGE+BMI+MEN+CURR_SMOKE+PREVAL_DIAB+BL_USE_RX_L+SYSTM+BP_TREAT + PC3,
                            data=meta(east_pseq),
                            x=TRUE),
                  pred = "PC3",
                  data=meta(east_pseq),
                  title = FALSE) +
  labs(x="", y="", subtitle="East", title = "") +
  coord_cartesian(ylim = c(-1.5, 2)) +
  theme_classic(20)


west_p <- plot_HR(fit=coxph(Surv(meta(west_pseq)$DEATH_AGEDIFF, meta(west_pseq)$DEATH) ~
                              BL_AGE+BMI+MEN+CURR_SMOKE+PREVAL_DIAB+BL_USE_RX_L+SYSTM+BP_TREAT + PC3,
                            data=meta(west_pseq),
                            x=TRUE),
                  pred = "PC3",
                  data=meta(west_pseq),
                  title = FALSE) +
  labs(x="", y="", subtitle="West", title = "") +
  coord_cartesian(ylim = c(-1.5, 2)) +
  scale_y_continuous(labels = NULL) +
  theme_classic(20)



p <- plot_grid(east_p, west_p, ncol = 2) +
  draw_label("PC3", x=0.5, y=  0, vjust=-0.5, angle= 0, size = 20) +
  draw_label("HR for Death", x=0.025, y=0.5,  angle= 90, size = 20)


# Save
if (!dir.exists(paste0(pltdir, "/figure_east_west_mortality/"))) {
  dir.create(paste0(pltdir, "/figure_east_west_mortality/"), recursive=TRUE)
}
png(filename = paste0(pltdir, "/figure_east_west_mortality/mortality_east_west_PC3.png"), 
    res = 300, width = 9, height = 4, units = "in")
p
dev.off()




