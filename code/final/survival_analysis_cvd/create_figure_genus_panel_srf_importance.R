## Data
pseq <- readRDS(file = genus_level_phyloseq_path)

# Add PCA components and diversity measures
pseq <- add_clr_abundances_to_meta_data(pseq)





## Genus panel 
genus_plots <- list()

# Yokenella
genus_plots[["Yokenella"]] <- plot_HR(fit=coxph(Surv(meta(pseq)$DEATH_AGEDIFF, meta(pseq)$DEATH) ~ 
                                                  BL_AGE+BMI+MEN+CURR_SMOKE+PREVAL_DIAB+BL_USE_RX_L + pspline(Yokenella_Bacteria),
                                                data=meta(pseq),
                                                x=TRUE),
                                      pred = "Yokenella_Bacteria",
                                      data=meta(pseq)) +
  labs(x="", y="", subtitle="Yokenella", title = "") +
  theme_classic(20) +
  coord_cartesian(ylim = c(-2.5, 4)) +
  scale_y_continuous(labels = NULL) 


# Butyrivibrio 
genus_plots[["Butyrivibrio"]] <- plot_HR(fit=coxph(Surv(meta(pseq)$DEATH_AGEDIFF, meta(pseq)$DEATH) ~ 
                                                     BL_AGE+BMI+MEN+CURR_SMOKE+PREVAL_DIAB+BL_USE_RX_L+SYSTM+BP_TREAT + Butyrivibrio_Bacteria,
                                                   data=meta(pseq),
                                                   x=TRUE),
                                         pred = "Butyrivibrio_Bacteria", data=meta(pseq)) +
  labs(x="", y="", subtitle="Butyrivibrio", title = "") +
  theme_classic(20) +
  coord_cartesian(ylim = c(-2.5, 4)) +
  # scale_y_continuous(labels = NULL) +
  scale_x_continuous(breaks = c(6, 9, 12, 15))

# Faecalitalea
genus_plots[["Faecalitalea"]] <- plot_HR(fit=coxph(Surv(meta(pseq)$DEATH_AGEDIFF, meta(pseq)$DEATH) ~ 
                                                     BL_AGE+BMI+MEN+CURR_SMOKE+PREVAL_DIAB+BL_USE_RX_L+SYSTM+BP_TREAT + Faecalitalea_Bacteria,
                                                   data=meta(pseq),
                                                   x=TRUE),
                                         pred = "Faecalitalea_Bacteria", data=meta(pseq)) +
  labs(x="", y="", subtitle="Faecalitalea", title = "") +
  theme_classic(20) +
  coord_cartesian(ylim = c(-2.5, 4)) +
  scale_y_continuous(labels = NULL)


# Fournierella
genus_plots[["Fournierella"]] <- plot_HR(fit=coxph(Surv(meta(pseq)$DEATH_AGEDIFF, meta(pseq)$DEATH) ~ BL_AGE+BMI+MEN+CURR_SMOKE+PREVAL_DIAB+BL_USE_RX_L + pspline(Fournierella_Bacteria),
                                                   data=meta(pseq), x=TRUE),
                                         pred = "Fournierella_Bacteria", data=meta(pseq)) +
  labs(x="", y="", subtitle="Fournierella", title = "") +
  theme_classic(22.5) +
  coord_cartesian(ylim = c(-2.5, 5)) +
  scale_y_continuous(labels = NULL)


genus_panel<- plot_grid(genus_plots[["Butyrivibrio"]],
                        genus_plots[["Faecalitalea"]], 
                        genus_plots[["Fournierella"]], 
                        genus_plots[["Yokenella"]], 
                        nrow = 1,
                        rel_widths = c(1.2, 1, 1, 1)) +
  draw_label("Abundance (CLR)", x=0.5, y=  0, vjust=-0.5, angle= 0, size = 20) +
  draw_label("HR for Death", x=0.025, y=0.5,  angle= 90, size = 20)





## SRF panel

rf_all <- readRDS(file = "output/raw_output/survival_random_forest_results.Rds")


# Data
importance_plot_data <- rf_all[["core"]]$importance %>% 
  m_neat(colnames = c("predictor", "importance")) %>%
  mutate(predictor = clean_genus_names(predictor, kingdom = FALSE), dummy = 1) %>% 
  arrange(desc(importance))


importance_panel <- importance_plot_data %>%
  head(n=20) %>%
  ggplot(aes(x = factor(predictor, levels = (.$predictor)), y = importance)) +
  geom_col(color = "black")+
  labs(x="", y="Importance \n", title = "") +
  theme_classic(22.5)+
  theme(axis.text.x = element_text(angle=45, hjust = 1))






p <- plot_grid(genus_panel,
               importance_panel,
               ncol = 1,
               labels = c("A", "B"),
               rel_heights = c(1, 1.3),
               label_size = 20)


if (!dir.exists(paste0(pltdir, "figure_genus_hazard_importance/"))) {
  dir.create(paste0(pltdir, "figure_genus_hazard_importance/"), recursive=TRUE)
}
png(filename = paste0(pltdir,"/figure_genus_hazard_importance/genus_hazard_importance.png"), 
    res = 300, width = 13, height = 10, units = "in")
p
dev.off()




