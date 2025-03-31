## Survival random forest 

# Load phyloseq
pseq <- readRDS(file = genus_level_phyloseq_path)

# Add CLR abundances
pseq <- dirty_genus_names(pseq)
pseq <- add_clr_abundances_to_meta_data(pseq)

data <- meta(pseq)

# Get core taxa
core_taxa <- core(pseq %>%
                    transform("compositional"), detection = .1/100, prevalence = 1/100) %>%
  taxa_names()



if(!exists("covariates")) covariates <- c("BL_AGE","BMI", "MEN",
                                          "CURR_SMOKE", "PREVAL_DIAB", "BL_USE_RX_L", 
                                          "SYSTM", "BP_TREAT")


# List of predictors: core, covariates, both core and covariates
rf_predictor_sets <- list(core = core_taxa,
                          covariates = covariates,
                          core_and_covariates = c(core_taxa, covariates))




# Run, forest, run!
rf_all <- lapply(names(rf_predictor_sets), function(x) {
  print(x)
  
  rf <- rfsrc(Surv(DEATH_AGEDIFF, DEATH) ~ .,
              data = meta(pseq)[, c("DEATH_AGEDIFF", "DEATH", rf_predictor_sets[[x]])],
              importance = TRUE)
  
  return(rf)
}) %>% 
  set_names(names(rf_predictor_sets))



# Save raw results
saveRDS(rf_all,
        file = "output/raw_output/survival_random_forest_results.Rds")

rf_all <- readRDS("output/raw_output/survival_random_forest_results.Rds")

# Get c-statistics
all_harrell_c <- sapply(rf_all, function(x) 1 - x$err.rate[length(x$err.rate)])

all_harrell_c <- all_harrell_c %>% 
  as.data.frame() %>% 
  m_neat %>% 
  set_colnames(c("predictor", "p_value"))



# Combine importances for the predictor sets
rsf_table <- data.frame(feature = rf_predictor_sets[["core_and_covariates"]])

for(pred in names(rf_predictor_sets)) {
  
  df <- rf_all[[pred]]$importance %>%
    as.data.frame() %>% 
    m_neat(colnames = c("feature", pred)) %>% 
    arrange()
  
  rsf_table <- full_join(rsf_table, df, "feature")
  
}


# Arrange to descing order 
rsf_table <- rsf_table %>%
  arrange(desc(core_and_covariates))

# Round
rsf_table[-1] <- round(rsf_table[-1], 5)

for(i in 2:4) {
  rsf_table[, i] <- ifelse(rsf_table[, i] == 0, 
                           "<0.00001",
                           as.character(rsf_table[, i]))
  
}




# Clean feature and column names
rsf_table <- rsf_table %>%
  mutate(feature = clean_genus_names(feature, kingdom = FALSE)) %>% 
  set_colnames(c("Feature", "Microbiome", "Covariates", "Microbiome and Covariates"))

rsf_table[rsf_table$Feature == "BLAGE", "Feature"] <-  "BL_AGE"




# Write tables
write.xlsx(rsf_table,
           file = paste0(out_tabledir , "/random_survival_forest_importance.xlsx"),
           colWidths = "auto")

write_delim(rsf_table,
            path = paste0(out_tabledir, "/random_survival_forest_importance.tsv"),
            delim = "\t")

write.xlsx(all_harrell_c,
           file = paste0(out_tabledir , "/random_survival_forest_c_statistic.xlsx"),
           colWidths = "auto")

