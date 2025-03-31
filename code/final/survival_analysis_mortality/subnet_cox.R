# Subnet Cox regression

print("subnet_Cox")


# Load phyloseq
pseq <- readRDS(file = genus_level_phyloseq_path)

# Subset pseqs
east_pseq <- pseq %>% 
  subset_samples(EAST == "EAST")
west_pseq <- pseq %>% 
  subset_samples(EAST == "WEST")

# Add CLR abundances
pseq <- dirty_genus_names(pseq)
pseq <- add_clr_abundances_to_meta_data(pseq)
sample_data(pseq) <- cbind(meta(pseq),
                           subnet_abundances(pseq,
                                             subnet_components,
                                             clr = TRUE))



if(!exists("alpha_level"))        alpha_level <- 0.05
if(!exists("status"))             status <- "DEATH"
if(!exists("time_to_event"))      time_to_event <- "DEATH_AGEDIFF"
if(!exists("splines"))            splines <- TRUE
if(!exists("normalize"))          normalize <- TRUE
if(!exists("test_ph_assumption")) test_ph_assumption <- FALSE

if(!exists("covariates")) covariates <- c("BL_AGE","BMI", "MEN",
                                          "CURR_SMOKE", "PREVAL_DIAB", "BL_USE_RX_L", 
                                          "SYSTM", "BP_TREAT")


subnet_abundances_clr <- meta(pseq)[, grep("subnet", colnames(meta(pseq)), value = TRUE)]

# predictors
subnet_predictors <- paste0("subnet_", 1:4)

# Data for the Cox wrapper
subnet_data <- subnet_abundances_clr %>% 
  as.data.frame() %>% 
  cbind(meta(pseq)[, c(covariates, "DEATH", "DEATH_AGEDIFF")])



# Cox regression
subnet_cox <- cox_wrapper(data = subnet_data,
                          predictors = subnet_predictors,
                          covariates = covariates,
                          alpha_level = alpha_level, ## FDR set to 1 to get all results to neat table
                          status = status,
                          time_to_event = time_to_event,
                          splines = splines,
                          normalize = normalize,
                          test_ph_assumption = FALSE)



output_table <- subnet_cox$neat_results %>% 
  mutate(`P (adjusted)` = round(`P (adjusted)`, 5))


# Write table
write.xlsx(output_table,
           file = paste0(out_tabledir, "cox_subnet.xlsx"), 
           colWidths = "auto")


## Write tsv table
write.table(output_table,
            file = paste0(out_tabledir, "cox_subnet.tsv"),
            sep = "\t", 
            row.names = FALSE)





## East vs West

# Add CLR abundances
east_pseq <- dirty_genus_names(east_pseq)
east_pseq <- add_clr_abundances_to_meta_data(east_pseq)
east_pseq <- cbind(meta(east_pseq), subnet_abundances(east_pseq, clr = TRUE))

west_pseq <- dirty_genus_names(west_pseq)
west_pseq <- add_clr_abundances_to_meta_data(west_pseq)
west_pseq <- cbind(meta(west_pseq), subnet_abundances(west_pseq, clr = TRUE))


# EAST

subnet_east_abundances_clr <- meta(east_pseq)[, grep("subnet", colnames(meta(east_pseq)), value = TRUE)]

# Data for the Cox wrapper
subnet_east_data <- subnet_east_abundances_clr %>% 
  as.data.frame() %>% 
  cbind(meta(east_pseq)[, c(covariates, "DEATH", "DEATH_AGEDIFF")])

  # Cox regression
subnet_east_cox <- cox_wrapper(data = subnet_east_data,
                          predictors = subnet_predictors,
                          covariates = covariates,
                          alpha_level = alpha_level,
                          status = status,
                          time_to_event = time_to_event,
                          splines = splines,
                          normalize = normalize,
                          test_ph_assumption = FALSE)

# WEST

subnet_west_abundances_clr <- meta(west_pseq)[, grep("subnet", colnames(meta(west_pseq)), value = TRUE)]

# Data for the Cox wrapper
subnet_west_data <- subnet_west_abundances_clr %>% 
  as.data.frame() %>% 
  cbind(meta(west_pseq)[, c(covariates, "DEATH", "DEATH_AGEDIFF")])

  # Cox regression
subnet_west_cox <- cox_wrapper(data = subnet_west_data,
                          predictors = subnet_predictors,
                          covariates = covariates,
                          alpha_level = alpha_level,
                          status = status,
                          time_to_event = time_to_event,
                          splines = splines,
                          normalize = normalize,
                          test_ph_assumption = FALSE)