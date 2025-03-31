## Survival random forest 

# Load phyloseq
pseq <- readRDS(file = genus_level_phyloseq_path)

# Add CLR abundances
pseq <- dirty_genus_names(pseq)
pseq <- add_clr_abundances_to_meta_data(pseq)

# Add causes of death
pseq <- get_causes_of_death(pseq)

data <- meta(pseq)
EAST_data <- data %>% filter(EAST == "EAST")
WEST_data <- data %>% filter(EAST == "WEST")

# Get core taxa IN THE EAST
core_taxa <- core(pseq %>%
                    subset_samples(EAST == "EAST") %>% 
                    transform("compositional"), detection = .1/100, prevalence = 1/100) %>%
  taxa_names()


if(!exists("covariates")) covariates <- c("BL_AGE","BMI", "MEN",
                                          "CURR_SMOKE", "PREVAL_DIAB", "BL_USE_RX_L", 
                                          "SYSTM", "BP_TREAT")


# List of predictors: core, covariates, both core and covariates
rf_predictor_sets <- list(core = core_taxa,
                          covariates = covariates,
                          core_and_covariates = c(core_taxa, covariates))



# Causes of death
causes <- as.character(unique(data$cause_of_death))
causes <- causes[!is.na(causes)]

# Cause specific Cox regression
cause_survival_forests <- lapply(c("All", causes), function(x) {
  
  # Create cause specific endpoint
  mydata <- data %>%
    mutate(CAUSE_DEATH = ifelse(cause_of_death != x,
                                0, 1)) %>%
    mutate(CAUSE_DEATH = ifelse(DEATH == 0, 0, CAUSE_DEATH))
  
  mydata_EAST <- mydata %>% filter(EAST == "EAST")
  mydata_WEST <- mydata %>% filter(EAST == "WEST")
  
  
  
  # Use differential status for causes
  my_status <- ifelse(x == "All", "DEATH", "CAUSE_DEATH")
  
  
  # Train survival random forest in EAST
  rf_train_EAST <- lapply(names(rf_predictor_sets), function(y) {
    print(y)
    
    srf_formula <- paste0("Surv(DEATH_AGEDIFF, ", my_status, ") ~ .") %>%
      as.formula
    
    rf <- rfsrc(srf_formula,
                data = mydata_EAST[, c("DEATH_AGEDIFF", my_status, rf_predictor_sets[[y]])],
                importance = TRUE)
    
    return(rf)
  }) %>% set_names(names(rf_predictor_sets))
  # Predict in the WEST
  rf_predict_WEST <- lapply(names(rf_predictor_sets), function(y) {
    print(y)
    
    rf <- predict(rf_train_EAST[[y]],
                  data = mydata_WEST[, c("DEATH_AGEDIFF", my_status, rf_predictor_sets[[y]])],
                  importance = TRUE)
    
    return(rf)
  }) %>% set_names(names(rf_predictor_sets))
  
  
  
  # Get prediction errors
  error_res <- lapply(rf_predict_WEST, function(rf) {
    
    last(rf$err.rate)
    
  }) %>%
    unlist %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    mutate(cause = x) %>%
    set_colnames(c("Predictor", "Error rate", "Cause of Death"))
  
  return(error_res)
  
}) %>%
  do.call(rbind, .)


output_table <- cause_survival_forests %>%  
  spread(`Predictor`, `Error rate`) %>% 
  set_colnames(c("Cause of Death", "Microbiome", "Covariates and Microbiome", "Covariates")) %>% 
  select(1, 2, 4, 3)

# Plot
cause_survival_forests %>% 
  ggplot(aes(x = `Cause of Death`, y = `Error rate`, color = Predictor)) + 
  geom_point()




# Write table
write.xlsx(output_table,
           file = paste0("output/tables/random_forest_EAST_WEST_accuracy.xlsx"),
           colWidths = "auto")


## Write tsv table
write.table(output_table,
            file = paste0("output/tables/random_forest_EAST_WEST_accuracy.tsv"),
            sep = "\t",
            row.names = FALSE)



# # Write table
# write.xlsx(output_table,
#            file = paste0("output/tables/random_forest_EAST_WEST_accuracy_norxj01.xlsx"),
#            colWidths = "auto")
# 
# 
# ## Write tsv table
# write.table(output_table,
#             file = paste0("output/tables/random_forest_EAST_WEST_accuracy_norxj01.tsv"),
#             sep = "\t",
#             row.names = FALSE)








