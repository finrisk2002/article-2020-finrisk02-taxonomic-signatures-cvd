source("code/final/R_functions.R")

run_cox_wrapper_with_covariates <- function(pseq.cox, 
  covariates = c("BL_AGE","BMI", "MEN", "CURR_SMOKE", "PREVAL_DIAB", "SYSTM", "BL_USE_RX_L", "BP_TREAT"),
  alpha_level =0.05 )
{
  taxa_names(pseq.cox) <- taxa_names_underscored(pseq.cox) %>% gsub(" ", "", .)



  cox_data <- create_cox_data(pseq.cox, covariates)


  cat("pre-cox printouts, n\n")
  cat(nsamples(pseq.cox))
  cat("\n")
  cat("complete cases\n")
  cat(sum(complete.cases(meta(pseq.cox)[,c(covariates, "DEATH")])))
  cat("\n")
  cat("number of deaths\n")
  cat(sum(meta(pseq.cox)$DEATH, na.rm=TRUE))
  cat("\n")

  cat("cox_data\n")
  cat(class(cox_data))
  cat("\n")
  cat(dim(cox_data))
  cat("\n")

  # cox_wrapper available at survival/survival_functions.R. TODO: check that no code duplicated
  pseq_cox_res <- cox_wrapper(data = cox_data,
                              predictors = taxa_names(pseq.cox),
                              covariates = covariates,
                              status = "DEATH",
                              time_to_event = "DEATH_AGEDIFF",
                              splines = TRUE,
                              normalize = TRUE,
                              test_ph_assumption = FALSE,
                              alpha_level = alpha_level)

  pseq_cox_res
}