# Cox regression: cause of death ~ Enterobacteriaceae

## Settings ****************** ####

print("cause_of_death_cox")
library(forestplot)
# Load phyloseq
pseq <- readRDS(file = genus_level_phyloseq_path)

# Edit genus names for a smoother R experience
pseq <- dirty_genus_names(pseq)

# Causes of death into pseq
pseq <- get_causes_of_death(pseq)

# if K_TPKS is not available, try
# pseq <- get_causes_of_death_ucod(pseq)

aggregate_taxa_clr_transform <- function(pseq,
                                         level = "Family") {
  
  
  tax_df <- tax_table(pseq) %>%
    as.data.frame() %>% 
    rownames_to_column()
  
  fams <- unique(tax_df  %>% pull(Family)) %>% 
    as.character
  
  
  aggregates <- lapply(fams, function(f) {
    
    genera <- tax_df %>% 
      filter(Family == f) %>% 
      pull(rowname) %>% 
      as.character()
    
    aggregated <- t(abundances(pseq))[, genera] %>% 
      as.data.frame %>% 
      rowSums()
    
    return(aggregated)
  }) %>%
    do.call(cbind, .)
  
  fams[is.na(fams)] <- "Other"
  colnames(aggregates) <- fams
  
  
  # CLR
  aggregated_clr <- aggregates %>% 
    t %>% 
    microbiome::transform("compositional") %>% 
    microbiome::transform("clr") %>% 
    t
  
  
  return(aggregated_clr)
  
}
sample_data(pseq) <- cbind(meta(pseq), aggregate_taxa_clr_transform(pseq, "Family"))


# Cox settings
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



causes <- as.character(unique(data$cause_of_death))
causes <- causes[!is.na(causes)]

# Cause specific Cox regression
cause_cox <- lapply(c("All", causes), function(x) {
  
  mydata <- data %>% 
    mutate(CAUSE_DEATH = ifelse(cause_of_death != x,
                                0, 1)) %>% 
    mutate(CAUSE_DEATH = ifelse(DEATH == 0, 0, CAUSE_DEATH))
  
  # Use differential status for causes
  my_status <- ifelse(x == "All", "DEATH", "CAUSE_DEATH")
  
  cox <- cox_wrapper_with_covariates(data = mydata,
                                     predictors = "Enterobacteriaceae",
                                     covariates = covariates,
                                     alpha_level = 1, ## here alpha = 1 to get non-significant results in neat format
                                     status = my_status,
                                     time_to_event = time_to_event,
                                     splines = splines,
                                     normalize = normalize)
  
  cox <- lapply(cox, function(df) {
    
    df %>% 
      mutate(Cause_Of_Death = x)
    
  })
  
  
  return(cox)
})


# Combine tables
cause_cox <- lapply(names(cause_cox[[1]]), function(tbl) {
  
  lapply(cause_cox, function(X) X[[tbl]]) %>% 
    do.call(rbind, .)
  
}) %>%
  set_names(names(cause_cox[[1]]))


## Table ********************* ####

# Get linear results in sweet format
output_table <- cause_cox$neat_results %>% 
  filter(Association == "linear", 
         Variable == "Enterobacteriaceae") %>% 
  mutate(FDR = p.adjust(P, "BH")) %>% 
  select(1:HR, FDR, everything(), -c(P))
# select(Predictor, HR, Cause_Of_Death, "P (adjusted)") %>% 
# mutate(HR = gsub("95% CI, ", "", HR)) %>% 
# spread(Cause_Of_Death, HR) %>% 
# mutate(Predictor = Predictor %>% 
#          gsub("_", " ", .) %>% 
#          gsub(" Bacteria", "", .))


# Write table
write.xlsx(output_table %>% 
             mutate(FDR = round(FDR, 5)),
           file = paste0("output/tables/cause_of_death_enterobacteriaceae.xlsx"), 
           colWidths = "auto")


## Write tsv table
write.table(output_table %>% 
              mutate(FDR = round(FDR, 5)),
            file = paste0("output/tables/cause_of_death_enterobacteriaceae.tsv"),
            sep = "\t", 
            row.names = FALSE)


## Forest plot *************** ####

## Forest Plot
plot_df <- cause_cox$results %>%
  filter(linearity == "linear", 
         variable == "Enterobacteriaceae") %>% 
  mutate(FDR = p.adjust(p, "BH"), 
         HR = paste0(round(exp(coef), 2), 
                     " (", round(exp(coef - 1.96*se_coef), 2), 
                     "-", round(exp(coef + 1.96*se_coef), 2), ")")) %>% 
  select(HR, FDR, Cause_Of_Death) %>% 
  mutate(HR_full = HR %>% 
           gsub("95% CI, ", "", .), 
         HR = gsub("95% CI, ", "", HR)) %>%
  separate(HR,
           into = c("HR", "CI"),
           sep = " ") %>% 
  separate(CI, 
           into = c("lower", "upper"), 
           sep = "-") %>%
  mutate(lower = gsub("\\(", "", lower) %>%
           as.numeric, 
         upper = gsub("\\)", "", upper) %>%
           as.numeric, 
         HR = HR %>% 
           as.numeric) %>% 
  mutate(n = ifelse(Cause_Of_Death == "Cancer", c(table(data$cause_of_death))["Cancer"],
                    ifelse(Cause_Of_Death == "Cardiovascular", c(table(data$cause_of_death))["Cardiovascular"],
                           ifelse(Cause_Of_Death == "Gastrointestinal", c(table(data$cause_of_death))["Gastrointestinal"],
                                  ifelse(Cause_Of_Death == "Neurological", c(table(data$cause_of_death))["Neurological"],
                                         ifelse(Cause_Of_Death == "Respiratory", c(table(data$cause_of_death))["Respiratory"],
                                                ifelse(Cause_Of_Death == "Other", c(table(data$cause_of_death))["Other"],
                                                       sum(table(data$cause_of_death))))))))) %>% 
  arrange(desc(HR))




numbers_df <- plot_df %>%  
  select(mean = HR, lower, upper) %>%
  rbind(rep(NA, 3), .)


text_df <- plot_df %>% 
  select(Cause_Of_Death, n, HR_full, FDR) %>% 
  mutate(FDR = round(FDR, 3))

# text_df <- rbind(c("Cause of Death", "Deaths", "HR (95% CI)", "FDR"),
#                  text_df)

## Manually add zeroes and use exponential notation to decimals for better visuals
text_df[text_df$Cause_Of_Death == "Gastrointestinal", "HR_full"] <- text_df[text_df$Cause_Of_Death == "Gastrointestinal", "HR_full"] %>%
  gsub("8 ", "80 ", .)

text_df[text_df$Cause_Of_Death == "Neurological", "HR_full"] <- text_df[text_df$Cause_Of_Death == "Neurological", "HR_full"] %>%
  gsub("9-", "90-", .)

text_df[text_df$Cause_Of_Death == "Other", "HR_full"] <- text_df[text_df$Cause_Of_Death == "Other", "HR_full"] %>%
  gsub("1.1", "1.10", .)


table_text <- list(
  # c(text_df$Cause_Of_Death),
  # c(text_df$n), 
  c("Cause of Death", text_df$Cause_Of_Death), 
  c("Deaths", text_df$n),
  c("HR", text_df$HR_full),
  c("FDR",
    append(append(list(expression("2.2x10"^-4)),
                  sprintf("%.3f", text_df$FDR[2:3])), 
           append(list(expression("9.5x10"^-5)),
                  sprintf("%.3f", text_df$FDR[5:7])))
  )
)



# png("forestplot.png", width=10, height=10, units = "in", res = 300)

#### Automatic saving doesn't work for whatever reason so do it manually
forestplot(table_text, 
           numbers_df,
           new_page = TRUE,
           is.summary = c(TRUE, rep(FALSE,7)),
           xlab = "Hazard Ratio",
           xticks = c(.75, 1, 2.5),
           x.ticks.digits = 0,
           xlog = TRUE,
           txt_gp = fpTxtGp(ticks  = gpar(cex = 1.5), 
                            xlab  = gpar(fontfamily = "", cex = 1.5),
                            label = gpar(cex = 1.5)),
           col = fpColors(box = "black",
                          line = "black",
                          summary = "royalblue"))




# print(p)
# 
# # # save plot
# dev.copy(png, "forestplot.png")
# dev.off()
# 
# 



# Cox regression: cause of death ~ Enterobacteriaceae

print("cause_of_death_cox")

pseq <- readRDS("input/data_work/phfinrisk_genus_all_drop50k_2018-12-21_nop_norxj01.rds")

# Load phyloseq
# pseq <- readRDS(file = genus_level_phyloseq_path)


# Function to get the most common causes of death, ICD class letter
get_causes_of_death <- function(x) {
  
  if(class(x) == "phyloseq") {
    causes <- meta(pseq)[, "K_TPKS"]
  } else {
    causes <- x
  }
  
  # Remove number part
  causes <- gsub("[0-9]+", "", causes) 
  
  for(i in 1:length(causes)) {
    if(is.na(causes[i])) {
      causes[i] <- NA
    } else if(causes[i] == "C") {
      causes[i] <- "Cancer"
    } else if(causes[i] == "G") {
      causes[i] <- "Neurological"
    } else if(causes[i] == "I") {
      causes[i] <- "Cardiovascular"
    } else if(causes[i] == "K") {
      causes[i] <- "Gastrointestinal"
    } else if(causes[i] == "J") {
      causes[i] <- "Respiratory"
    } else if(!is.na(causes[i])) {
      causes[i] <- "Other"
    }
  }
  
  if(class(x) == "phyloseq") {
    
    cause_of_death <- causes
    
    sample_data(x) <- cbind(meta(x), cause_of_death)
    return(x)
    
  } else {
    return(causes)
  }
  
  
}

# Add causes to pseq
pseq <- get_causes_of_death(pseq)

# Edit genus names for a smoother R experience
pseq <- dirty_genus_names(pseq)


# Aggregate_taxa from microbiome doesn't work for whatever reason
aggregate_taxa_clr_transform <- function(pseq,
                                       level = "Family") {
  
  tax_df <- tax_table(pseq) %>% as.data.frame()
  otu_df <- t(abundances(pseq))
  
  # Get genera in the given group
  levels <- tax_df[, level] %>% 
    as.character() %>% 
    unique() 
  
  # Genera in the specific group
  aggregates <- lapply(levels, function(L) {
    # print(L)
    genera <- tax_df[tax_df[, level] == L, ] %>% 
      drop_na %>% 
      rownames() %>% 
      as.character()
    
    aggregate <- otu_df[, genera] %>% 
      as.data.frame() %>% 
      rowSums()
    
    return(aggregate)
    
  }) %>%
    do.call(cbind, .)
    
  colnames(aggregates) <- replace(levels, is.na(levels), "Other")

  # CLR transform
  otu_df_clr <- aggregates %>% 
    t %>% 
    microbiome::transform("compositional") %>% 
    microbiome::transform("clr") %>% 
    t
  
  return(otu_df_clr)
}


sample_data(pseq) <- cbind(meta(pseq), aggregate_taxa_clr_transform(pseq, "Family"))


# Cox settings
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



causes <- as.character(unique(data$cause_of_death))
causes <- causes[!is.na(causes)]

# Cause specific Cox regression
cause_cox <- lapply(causes, function(x) {
  
  mydata <- data %>% 
    filter(cause_of_death == x |
             is.na(cause_of_death))
  
  cox <- cox_wrapper(data = mydata,
                           predictors = "Enterobacteriaceae",
                           covariates = covariates,
                           alpha_level = 1, ## here alpha = 1 to get non-significant results in neat format
                           status = status,
                           time_to_event = time_to_event,
                           splines = splines,
                           normalize = normalize,
                           test_ph_assumption = FALSE)
  
  cox <- lapply(cox, function(df) {
    
    df %>% 
      mutate(Cause_Of_Death = x)
    
  })
  
  
  return(cox)
})


# Combine tables
cause_cox <- lapply(names(cause_cox[[1]]), function(tbl) {
  
  lapply(cause_cox, function(X) X[[tbl]]) %>% 
    do.call(rbind, .)
  
}) %>%
  set_names(names(cause_cox[[1]]))


# Get linear results in sweet format
output_table <- cause_cox$neat_results %>% 
  filter(Association == "linear") 
  # select(Predictor, HR, Cause_Of_Death, "P (adjusted)") %>% 
  # mutate(HR = gsub("95% CI, ", "", HR)) %>% 
  # spread(Cause_Of_Death, HR) %>% 
  # mutate(Predictor = Predictor %>% 
  #          gsub("_", " ", .) %>% 
  #          gsub(" Bacteria", "", .))


# Write table
write.xlsx(output_table,
           file = paste0("output/tables/.xlsx"), 
           colWidths = "auto")


## Write tsv table
write.table(output_table,
            file = paste0("output/tables/.tsv"),
            sep = "\t", 
            row.names = FALSE)


## Forest Plot

plot_df <- output_table %>% 
  select(HR, "P (adjusted)", Cause_Of_Death) %>% 
  mutate(HR = gsub("95% CI, ", "", HR)) %>%
  separate(HR,
           into = c("HR", "CI"),
           sep = " ") %>% 
  separate(CI, 
           into = c("lower", "upper"), 
           sep = "-") %>%
  mutate(lower = gsub("\\(", "", lower) %>%
           as.numeric, 
         upper = gsub("\\)", "", upper) %>%
           as.numeric, 
         HR = HR %>% 
           as.numeric) %>% 
  mutate(n = ifelse(Cause_Of_Death == "Cancer", c(table(data$cause_of_death))["Cancer"],
                    ifelse(Cause_Of_Death == "Cardiovascular", c(table(data$cause_of_death))["Cardiovascular"],
                           ifelse(Cause_Of_Death == "Gastrointestinal", c(table(data$cause_of_death))["Gastrointestinal"],
                                  ifelse(Cause_Of_Death == "Neurological", c(table(data$cause_of_death))["Neurological"],
                                         ifelse(Cause_Of_Death == "Respiratory", c(table(data$cause_of_death))["Respiratory"],
                                                c(table(data$cause_of_death))["Other"])))))) %>% 
  arrange(desc(HR))




numbers_df <- plot_df %>% 
  select(mean = HR, lower, upper) %>% 
  rbind(rep(NA, 3), .)


text_df <- plot_df %>% 
  select(Cause_Of_Death, n, HR, P_value = "P (adjusted)") %>% 
  mutate(P_value = round(P_value, 5))

text_df <- rbind(c("Cause", "Deaths", "HR", "P-value"),
                 text_df)


# png("forestplot.png", width=10, height=10, units = "in", res = 300)

forestplot(text_df, 
           numbers_df,
           new_page = TRUE,
           is.summary = c(TRUE, rep(FALSE,6)),
           xticks.digits = 1, 
           xlog = TRUE, 
           col = fpColors(box = "black",
                          line = "black",
                          summary = "royalblue"))

# # save plot
# dev.copy(png, "forestplot.png")
# dev.off()
# 
# 



