
print("PC3_checks_for_review")


# Load phyloseq
pseq <- readRDS(file = genus_level_phyloseq_path)
species_pseq <- readRDS(species_level_phyloseq_path)

# Add PCA components and diversity measures
pseq <- add_pca_components_to_meta_data(pseq, species_pseq, n_axes = 3)

pseq <- add_diversities_to_meta_data(pseq, species_pseq)

pseq <- get_causes_of_death(pseq)

data <- meta(pseq) 





if(!exists("alpha_level"))        alpha_level <- 0.05
if(!exists("status"))             status <- "DEATH"
if(!exists("time_to_event"))      time_to_event <- "DEATH_AGEDIFF"
if(!exists("splines"))            splines <- TRUE
if(!exists("normalize"))          normalize <- TRUE
if(!exists("test_ph_assumption")) test_ph_assumption <- FALSE

if(!exists("covariates")) covariates <- c("BL_AGE","BMI", "MEN",
                                          "CURR_SMOKE", "PREVAL_DIAB", "BL_USE_RX_L", 
                                          "SYSTM", "BP_TREAT")


PC3_pointer <- c("PC3")


var_description_df <- data.frame(Variable = c("BL_AGE",
                                              "MENMale",
                                              "CURR_SMOKE1", 
                                              "PREVAL_DIAB1", 
                                              "SYSTM", 
                                              "BL_USE_RX_L1", 
                                              "BMI", 
                                              "PC3", 
                                              "BP_TREAT", 
                                              "HFC_score"), 
                                 Description = c("Baseline Age", 
                                                 "Sex", 
                                                 "Smoking", 
                                                 "Diabetes", 
                                                 "Systolic Blood Pressure", 
                                                 "Antineoplastic or immunomodulating agents", 
                                                 "BMI", 
                                                 "PC3",
                                                 "Antihypertensive Medication",
                                                 "Healthy Diet Index"))





## PC3 Cox w covariates ******************* ####


# Cox regression
PC3_cox_full <- cox_wrapper_with_covariates(data = data,
                                            predictors = PC3_pointer,
                                            covariates = covariates,
                                            alpha_level = 1,   # Here FRD level set to 1 so the neat table contains non-significant results as requested
                                            status = status,
                                            time_to_event = time_to_event,
                                            splines = splines,
                                            normalize = normalize)




output_table <- PC3_cox_full$neat_results %>% 
  left_join(var_description_df, by = "Variable") %>% 
  mutate(P = ifelse(P < 1e-5, "<0.00001", as.character(round(P, 5)))) %>%
  select(Variable = Description, 
         Coefficient,
         "Coefficient SE", 
         HR, 
         P, 
         test_stat = "Test Statistic Value", 
         "Test Statistic", 
         Association) %>% 
  mutate(test_stat = round(test_stat, 3)) %>% 
  set_colnames(replace(colnames(.), 7, "Test Statistic"))


# Write table
write.xlsx(output_table,
           file = paste0("output/tables/PC3_cox_all_covariates.xlsx"), 
           colWidths = "auto")


## Write tsv table
write.table(output_table,
            file = paste0("output/tables/PC3_cox_all_covariates.tsv"),
            sep = "\t", 
            row.names = FALSE)


# # Write table
# write.xlsx(output_table,
#            file = paste0("output/tables/PC3_cox_all_covariates_norxj01.xlsx"), 
#            colWidths = "auto")
# 
# 
# ## Write tsv table
# write.table(output_table,
#             file = paste0("output/tables/PC3_cox_all_covariates_norxj01.tsv"),
#             sep = "\t", 
#             row.names = FALSE)



## PC3 Cox Diet index ********************* ####


# Cox regression
PC3_cox_full_HFC_score <- cox_wrapper_with_covariates(data = data,
                                                      predictors = PC3_pointer,
                                                      covariates = c(covariates, "HFC_score"),
                                                      alpha_level = 1,   # Here FRD level set to 1 so the neat table contains non-significant results as requested
                                                      status = status,
                                                      time_to_event = time_to_event,
                                                      splines = splines,
                                                      normalize = normalize)



output_table_HFC <- PC3_cox_full_HFC_score$neat_results %>% 
  left_join(var_description_df, by = "Variable") %>% 
  select(Variable, 
         Description, 
         Coefficient,
         "Coefficient SE", 
         HR, 
         P, 
         "Test Statistic Value", 
         "Test Statistic", 
         Association)



# Write table
write.xlsx(output_table_HFC,
           file = paste0("output/tables/PC3_cox_all_covariates_HFC_score.xlsx"), 
           colWidths = "auto")


## Write tsv table
write.table(output_table_HFC,
            file = paste0("output/tables/PC3_cox_all_covariates_HFC_score.tsv"),
            sep = "\t", 
            row.names = FALSE)

# # Write table
# write.xlsx(output_table_HFC,
#            file = paste0("output/tables/PC3_cox_all_covariates_HFC_score_norxj01.xlsx"), 
#            colWidths = "auto")
# 
# 
# ## Write tsv table
# write.table(output_table_HFC,
#             file = paste0("output/tables/PC3_cox_all_covariates_HFC_score_norxj01.tsv"),
#             sep = "\t", 
#             row.names = FALSE)
## PC3 exclude prevalent cancer *********** ####

cancer_data <- data %>% 
  filter(PREVAL_CR_ANYCANC != 1)


PC3_cox_cancer <- cox_wrapper(data = cancer_data,
                              predictors = PC3_pointer,
                              covariates = covariates,
                              alpha_level = 1,   # Here FRD level set to 1 so the neat table contains non-significant results as requested
                              status = status,
                              time_to_event = time_to_event,
                              splines = splines,
                              normalize = normalize)


cancer_output_table <- PC3_cox_cancer$neat_results %>% 
  filter(Association == "linear")



# Write table
write.xlsx(cancer_output_table,
           file = paste0("output/tables/PC3_cox_cancer.xlsx"), 
           colWidths = "auto")


## Write tsv table
write.table(cancer_output_table,
            file = paste0("output/tables/PC3_cox_cancer.tsv"),
            sep = "\t", 
            row.names = FALSE)
## PC3 Cox ******************************** ####

if(!exists("alpha_level"))        alpha_level <- 0.05
if(!exists("status"))             status <- "DEATH"
if(!exists("time_to_event"))      time_to_event <- "DEATH_AGEDIFF"
if(!exists("splines"))            splines <- TRUE
if(!exists("normalize"))          normalize <- TRUE
if(!exists("test_ph_assumption")) test_ph_assumption <- FALSE

if(!exists("covariates")) covariates <- c("BL_AGE","BMI", "MEN",
                                          "CURR_SMOKE", "PREVAL_DIAB", "BL_USE_RX_L", 
                                          "SYSTM", "BP_TREAT")


PC3_pointer <- c("PC3")



# Cox regression
PC3_cox_full <- cox_wrapper_with_covariates(data = data,
                        predictors = PC3,
                        covariates = covariates,
                        alpha_level = 1,   # Here FRD level set to 1 so the neat table contains non-significant results as requested
                        status = status,
                        time_to_event = time_to_event,
                        splines = splines,
                        normalize = normalize)



output_table <- PC3_cox_full$neat_results %>% 
  mutate(Description = c("Baseline Age", 
                         "Sex", 
                         "Smoking", 
                         "Diabetes", 
                         "Systolic Blood Pressure", 
                         "Antineoplastic or immunomodulating agents", 
                         "BMI", 
                         "PC3",
                         "Antihypertensive Medication")) %>% 
  select(Variable, 
         Description, 
         Coefficient,
         "Coefficient SE", 
         HR, 
         P, 
         "Test Statistic Value", 
         "Test Statistic", 
         Association)


# Write table
write.xlsx(output_table,
           file = paste0("output/tables/PC3_cox_all_covariates.xlsx"), 
           colWidths = "auto")


## Write tsv table
write.table(output_table,
            file = paste0("output/tables/PC3_cox_all_covariates.tsv"),
            sep = "\t", 
            row.names = FALSE)






