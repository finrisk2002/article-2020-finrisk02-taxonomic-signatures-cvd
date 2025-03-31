pseq <- readRDS(file = genus_level_phyloseq_path)

if(!exists("alpha_level"))        alpha_level <- 0.05
if(!exists("status"))             status <- "DEATH"
if(!exists("time_to_event"))      time_to_event <- "DEATH_AGEDIFF"
if(!exists("splines"))            splines <- TRUE
if(!exists("normalize"))          normalize <- TRUE
if(!exists("test_ph_assumption")) test_ph_assumption <- FALSE

if(!exists("covariates")) covariates <- c("BL_AGE","BMI", "MEN",
                                          "CURR_SMOKE", "PREVAL_DIAB", "BL_USE_RX_L", 
                                          "SYSTM", "BP_TREAT")



## Load and edit data **************** ####
functional_A <- readRDS(functional_activities_path)



## metabolite names are difficult to handle so create aliases
metabolite_names <- list(modules = data.frame(Metabolite = functional_A[[1]] %>%
                                                rownames,
                                              Alias = paste0("m", 1:nrow(functional_A[[1]]))),
                         pathways = data.frame(Metabolite = functional_A[[2]] %>%
                                                 rownames,
                                               Alias = paste0("p", 1:nrow(functional_A[[2]]))),
                         ko = data.frame(Metabolite = functional_A[[3]] %>%
                                           rownames,
                                         Alias = paste0("k", 1:nrow(functional_A[[3]])))
)


# Set alias colnames
functional_A <- lapply(names(functional_A), function(x) {
  
  df <- functional_A[[x]]
  
  df <- df %>%
    t %>% 
    set_colnames(metabolite_names[[x]]$Alias)
  
  
  return(df)
  
  
}) %>% set_names(names(functional_A))


# Add nuggets and log10 transform
modules <- functional_A[["modules"]] %>%
  apply(2, as.numeric) %>%
  as.data.frame() %>%
  + min(.[.!= 0]) %>% 
  log10()

pathways <- functional_A[["pathways"]] %>%
  apply(2, as.numeric) %>% 
  as.data.frame() %>%
  + min(.[.!= 0]) %>% 
  log10()

ko <- functional_A[["ko"]] %>%
  apply(2, as.numeric) %>%
  as.data.frame() %>%
  + min(.[.!= 0]) %>%
  log10()




# Ignore variables from ko that have only 2 different values + some additional features 
# that otherwise would cause trouble downstream
ko_n <- lapply(colnames(ko), function(x) {
  
  print(x)
  
  ko[, x] %>%
    table %>%
    as.data.frame() %>%
    nrow
  
}) %>% unlist %>% set_names(colnames(ko))
ko <- ko[, !(colnames(ko) %in% c(names(ko_n[ko_n == 2]), "k6441"))]


# Data frames used in Cox regression
ko <- cbind(ko, meta(pseq)[, c(covariates, "DEATH", "DEATH_AGEDIFF")]) %>% as.data.frame()


## Cox ******************************* ####

predictors <-  colnames(ko)[1:(ncol(ko) - length(covariates) - 2)]

ko_cox <- cox_regression(data = ko,
                         predictors,
                         covariates,
                         status = "DEATH",
                         time_to_event = "DEATH_AGEDIFF",
                         splines,
                         normalize)

# Remove pathological cases
ko_remove <- c("k126", "k1676", "k126", "k3401", "k6782", "k6767", "k6680",
               "k6535", "k6482", "k6471", "k6441", "k6440", "k6375", "k6155",
               "k6141", "k6137", "k6136", "k6123", "k6122", "k5954","k5507",
               "k5416", "k5415", "k5130", "k5078", "k5129", "k4218", "k4193",
               "k3834", "k3680", "k3640", "k3627", "k3640", "k3622", "k3511",
               "k3640")

predictors <- predictors[!(predictors %in% ko_remove)]


ko_cox_results <- cox_results(data = ko,
                              predictors = predictors,
                              covariates = covariates,
                              status = "DEATH",
                              time_to_event = "DEATH_AGEDIFF",
                              alpha_level = alpha_level,
                              splines = splines,
                              fit_list = ko_cox)


ko_cox_results[["neat_results"]] <- ko_cox_results[["neat_results"]] %>%
  mutate(Alias = Predictor) %>% 
  inner_join(metabolite_names[["ko"]], by = "Alias") %>% 
  select(Metabolite, Alias, Coefficient, "Coefficient SE", HR, "P (adjusted)", "Test Statistic Value", "Test Statistic", "Association")

ko_cox_results[["results"]] <- ko_cox_results[["results"]] %>%
  mutate(Alias = predictor) %>% 
  inner_join(metabolite_names[["ko"]], by = "Alias") %>% 
  select(Metabolite, Alias, coef, se_coef, PH, p, p_adj, posneg, association)


ko_cox_results[["linear_results"]] <- ko_cox_results[["linear_results"]] %>%
  mutate(Alias = predictor) %>% 
  inner_join(metabolite_names[["ko"]], by = "Alias") %>% 
  select(Metabolite, Alias, coef, se_coef, PH, p, p_adj, posneg, association)


## Modify KO Group results table

# get ko group terms
ko_group_terms <- read_file(file = ko_metabolite_terms_path)

ko_group_terms <- ko_group_terms %>%
  str_split("\n") %>% 
  sapply(function(x) str_split(x, "   ")) %>%
  do.call(rbind, .) %>% 
  as_tibble() %>% 
  separate(V1, sep = ":", into = c("KO", "KO_term")) %>% 
  separate(V2, sep = ";", into = c("code", "KO_group")) %>% 
  select(-c(KO, code)) %>% 
  mutate(KO_group = trimws(KO_group))



# Change column class to character
ko_cox_results[["neat_results"]]$Alias <- ko_cox_results[["neat_results"]]$Metabolite %>% as.character()
ko_cox_results[["neat_results"]]$Metabolite <- ko_cox_results[["neat_results"]]$Metabolite %>% as.character()

# Add metabolite term to neat results
for(i in 1:nrow(ko_cox_results[["neat_results"]])) {
  
  met <- ko_cox_results[["neat_results"]][i, "Alias"] %>% unlist
  
  text <- ko_group_terms[ko_group_terms$KO_term == met, "KO_group"]
  
  ko_cox_results[["neat_results"]][i, "Metabolite"] <- ifelse(nrow(text) == 0, NA, text)
}




# Change column class to character
ko_cox_results[["linear_results"]]$Alias <- ko_cox_results[["linear_results"]]$Metabolite %>% as.character()
ko_cox_results[["linear_results"]]$Metabolite <- ko_cox_results[["linear_results"]]$Metabolite %>% as.character()

# Add metabolite term to linear results
for(i in 1:nrow(ko_cox_results[["linear_results"]])) {
  
  met <- ko_cox_results[["linear_results"]][i, "Alias"] %>% unlist
  
  text <- ko_group_terms[ko_group_terms$KO_term == met, "KO_group"]
  
  ko_cox_results[["linear_results"]][i, "Metabolite"] <- ifelse(nrow(text) == 0, NA, text)
}



## Write tables ********************** ####

write.xlsx(ko_cox_results[["neat_results"]]  %>%
             filter(Association == "linear"), 
           file = "output/tables/cox_functional_ko.xlsx")
write.table(ko_cox_results[["neat_results"]]  %>% 
              filter(Association == "linear"),
            file = "output/tables/cox_functional_ko.tsv", 
            sep = "\t",
            row.names = FALSE)

write.xlsx(ko_cox_results[["linear_results"]], file = "output/tables/cox_linear_functional_ko.xlsx")






