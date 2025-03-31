# Cox regression for PCA 1-3, Shannon diversity and Observed richness

print("PCA_cox")


# Load phyloseq
pseq <- readRDS(file = genus_level_phyloseq_path)
species_pseq <- readRDS(species_level_phyloseq_path)

# Add PCA components and diversity measures
pseq <- add_pca_components_to_meta_data(pseq, species_pseq, n_axes = 3)
pseq <- add_diversities_to_meta_data(pseq)

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


# Predictors
misc_predictors <- c("PC1", "PC2", "PC3", "Shannon", "Observed")

# Cox regression
misc_cox <- cox_wrapper(data = data,
                        predictors = misc_predictors,
                        covariates = covariates,
                        alpha_level = 1,   # Here FRD level set to 1 so the neat table contains non-significant results as requested
                        status = status,
                        time_to_event = time_to_event,
                        splines = splines,
                        normalize = normalize,
                        test_ph_assumption = FALSE)



output_table <- misc_cox$neat_results %>% 
  filter(Association == "linear")


# Write table
write.xlsx(output_table,
           file = paste0(out_tabledir , "/cox_PCA_diversity.xlsx"), 
           colWidths = "auto")


## Write tsv table
write.table(output_table,
            file = paste0(out_tabledir , "/cox_PCA_diversity.tsv"),
            sep = "\t", 
            row.names = FALSE)

