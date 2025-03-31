## Survival random forest cross validation

# Load phyloseq
pseq <- readRDS(file = genus_level_phyloseq_path)

# Add CLR abundances
pseq <- dirty_genus_names(pseq)
pseq <- add_clr_abundances_to_meta_data(pseq)

data <- meta(pseq)

# Get core taxa
core_taxa <- core(pseq %>%
                    transform("compositional"),
                  detection = .1/100,
                  prevalence = 1/100) %>%
  taxa_names()



if(!exists("covariates")) covariates <- c("BL_AGE","BMI", "MEN",
                                          "CURR_SMOKE", "PREVAL_DIAB", "BL_USE_RX_L", 
                                          "SYSTM", "BP_TREAT")


# List of predictors: core, covariates, both core and covariates
rf_predictor_sets <- list(core = core_taxa,
                          covariates = covariates,
                          core_and_covariates = c(core_taxa, covariates))




# Cross validation
# NOTA BENE! This takes a long time 
rf_cv <-  random_survival_forest_cv(data = meta(pseq),
                                    predictor_list = rf_predictor_sets,
                                    response_var = "DEATH",
                                    response_time = "DEATH_AGEDIFF",
                                    n_folds = 5,
                                    seed = seed)


# Save raw output
save(rf_cv,
     file = paste0(datadir, "survival_random_forest_cross_validation.Rdata"))



# Get c-values
cv_harrell_c <- lapply(names(rf_cv), function(pred) {
  
  sapply(1:5, function(i) {
    rf_cv[[pred]][["cv"]][[i]][["harrell_c"]]
  }) %>%
    as.data.frame() %>% 
    cbind(., 1:5, rep(pred, 5)) %>% 
    set_colnames(c("harrell_c", "fold", "predictor"))
  
}) %>% do.call(rbind, .)


# Paired t-test: difference between feature sets

# unique levels
lvls <- cv_harrell_c$predictor %>% unique
harrell_t_test <- c()
harrell_t_test[1] <- t.test(cv_harrell_c[cv_harrell_c$predictor == lvls[1], ]$harrell_c,
                            cv_harrell_c[cv_harrell_c$predictor == lvls[2], ]$harrell_c,
                            paired = TRUE)$p.value
harrell_t_test[2] <- t.test(cv_harrell_c[cv_harrell_c$predictor == lvls[1], ]$harrell_c,
                            cv_harrell_c[cv_harrell_c$predictor == lvls[3], ]$harrell_c,
                            paired = TRUE)$p.value
harrell_t_test[3] <- t.test(cv_harrell_c[cv_harrell_c$predictor == lvls[2], ]$harrell_c,
                            cv_harrell_c[cv_harrell_c$predictor == lvls[3], ]$harrell_c,
                            paired = TRUE)$p.value

names(harrell_t_test) <- c("core_covariates",
                           "core_core_and_covariates",
                           "covariates_core_and_covariates")

# Into data frame for easier saving and downstream usage
harrell_t_test <- harrell_t_test %>% 
  as.data.frame() %>% 
  m_neat %>% 
  set_colnames(c("predictor", "p_value"))


# Write tables
write.xlsx(cv_harrell_c, file = paste0(out_tabledir , "survival_random_forest_cross_validation_c_statistics.xlsx"))
write.xlsx(harrell_t_test, file = paste0(out_tabledir , "/c_statistic_t_test_p_values.xlsx"))
