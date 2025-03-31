# Genus level Cox regression 
print("genus_Cox")


# Load phyloseq
pseq <- readRDS(file = genus_level_phyloseq_path)

# Edit genus names for a smoother R experience
pseq <- dirty_genus_names(pseq)

# Add CLR abundances
pseq <- add_clr_abundances_to_meta_data(pseq)

# Data used in subsequent analyses
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


# Predictors: genus names
genus_predictors <- taxa_names(pseq)

genus_cox <- cox_wrapper(data = data,
                        predictors = genus_predictors,
                        covariates = covariates,
                        alpha_level = alpha_level,  
                        status = status,
                        time_to_event = time_to_event,
                        splines = splines,
                        normalize = normalize,
                        test_ph_assumption = FALSE)



# Edit results
output_table <-  genus_cox[["neat_results"]] %>% 
  mutate(Genus = Predictor) %>%
  select(one_of("Genus", "Coefficient", "Coefficient SE", "HR",
                "P (adjusted)", "Test Statistic Value", "Test Statistic", "Association"))


output_table <- output_table %>% 
  mutate(Predictor = clean_genus_names(Genus)) %>% 
  mutate(Genus = Predictor) %>% 
  select(-Predictor)


## Round P values
output_table$`P (adjusted)` <- round(output_table$`P (adjusted)`, 5)



## Write xlsx table
write.xlsx(output_table,
           file = paste0(out_tabledir , "cox_genus.xlsx"))

## Write tsv table
write.table(output_table,
           file = paste0(out_tabledir , "cox_genus.tsv"),
           sep = "\t", 
           row.names = FALSE)










