# Cox regression: cause of death ~ Enterobacteriaceae

## Settings ************************ ####

print("cause_of_death_cox")
library(forestplot)
# Load phyloseq
pseq <- readRDS(file = genus_level_phyloseq_path)

# Edit genus names for a smoother R experience
pseq <- dirty_genus_names(pseq)

# Causes of death
pseq <- get_causes_of_death(pseq)

# Aggregate to family level
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

# Add Enterobaacteriaceae classes: lowest/highest 25%
sample_data(pseq) <- cbind(meta(pseq),
                           Enterobacteriaceae_Quartile = ntile(meta(pseq)$Enterobacteriaceae, 4) %>% 
                             as.factor())

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

# Causes of death
causes <- as.character(unique(data$cause_of_death))
causes <- causes[!is.na(causes)]




## Cox ***************************** ####

# Cause specific Cox regression
cause_cox <- lapply(c("All", causes), function(x) {
  
  mydata <- data %>% 
    mutate(CAUSE_DEATH = ifelse(cause_of_death != x,
                                0, 1)) %>% 
    mutate(CAUSE_DEATH = ifelse(DEATH == 0, 0, CAUSE_DEATH))
  
  # Use differential status for causes
  my_status <- ifelse(x == "All", "DEATH", "CAUSE_DEATH")
  
  
  cox <- cox_wrapper_categorical(data = mydata,
                                     predictors = "Enterobacteriaceae_Quartile",
                                     covariates = covariates,
                                     alpha_level = 1,   # Here FRD level set to 1 so the neat table contains non-significant results as requested
                                     status = my_status,
                                     time_to_event = time_to_event,
                                     normalize = normalize)
  
  
  cox$results %>% 
    filter(linearity == "linear", variable == "Enterobacteriaceae_Quartile4") %>% 
    select(coef, se_coef, p) %>% 
    mutate(lower_2.5 = coef - 1.96*`se_coef`, 
           upper_97.5 = coef + 1.96*`se_coef`) %>% 
    mutate(cause = x)
  
  
}) %>% 
  do.call(rbind, .)


## Quantile binned forest plot ***** ####

## Forest Plot
plot_df <- cause_cox %>%
  mutate(FDR = p.adjust(p, "BH"), 
         HR = paste0(round(exp(coef), 2), 
                     " (", round(exp(coef - 1.96*se_coef), 2), 
                     "-", round(exp(coef + 1.96*se_coef), 2), ")")) %>% 
  select(HR, FDR, Cause_Of_Death = cause) %>% 
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



# Table for numbers
numbers_df <- plot_df %>%  
  select(mean = HR, lower, upper) %>%
  rbind(rep(NA, 3), .)

# Table for text
text_df <- plot_df %>% 
  select(Cause_Of_Death, n, HR_full, FDR) %>% 
  mutate(FDR = as.character(round(FDR, 3)))



## Manually add zeroes and use exponential notation to decimals for better visuals
text_df[text_df$Cause_Of_Death == "Gastrointestinal", "HR_full"] <- text_df[text_df$Cause_Of_Death == "Gastrointestinal", "HR_full"] %>%
  gsub("4 ", "4.00 ", .)

text_df[text_df$Cause_Of_Death == "Other", "HR_full"] <- text_df[text_df$Cause_Of_Death == "Other", "HR_full"] %>%
  gsub("1.2", "1.20", .)

text_df[text_df$Cause_Of_Death == "Other", "HR_full"] <- text_df[text_df$Cause_Of_Death == "Cancer", "HR_full"] %>%
  gsub("0.7", "0.70", .)

text_df[text_df$Cause_Of_Death == "Cancer", "FDR"] <- text_df[text_df$Cause_Of_Death == "Cancer", "FDR"] %>%
  gsub("0.04", "0.040", .)


table_text <- list(
  # c(text_df$Cause_Of_Death),
  # c(text_df$n), 
  c("Cause of Death", text_df$Cause_Of_Death), 
  c("Deaths", text_df$n),
  c("HR (95% CI)", text_df$HR_full),
  c("FDR", text_df$FDR)
  # c("FDR",
  #   append(append(list(expression("2.2x10"^-4)),
  #                 sprintf("%.3f", text_df$FDR[2:3])), 
  #          append(list(expression("9.5x10"^-5)),
  #                 sprintf("%.3f", text_df$FDR[5:7])))
  # )
)



# png("forestplot.png", width=10, height=10, units = "in", res = 300)

#### Automatic saving doesn't work for whatever reason so do it manually
forestplot(table_text, 
           numbers_df,
           new_page = TRUE,
           is.summary = c(TRUE, rep(FALSE,7)),
           xlab = "Hazard Ratio",
           xticks = c(0.25, .5, 1, 2, 4, 8, 16),
           x.ticks.digits = 0,
           xlog = TRUE,
           txt_gp = fpTxtGp(ticks  = gpar(cex = 1.5), 
                            xlab  = gpar(fontfamily = "", cex = 1.5),
                            label = gpar(cex = 1.5)),
           col = fpColors(box = "black",
                          line = "black",
                          summary = "royalblue"))

