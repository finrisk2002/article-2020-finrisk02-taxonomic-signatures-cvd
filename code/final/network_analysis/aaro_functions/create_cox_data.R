# helper function
create_cox_data <- function(pseq_tmp_cox, covariates) {

  pseq_tmp2_cox <- pseq_tmp_cox
  death_vars <- meta(pseq_tmp2_cox)[, c("DEATH", "DEATH_AGEDIFF")]

  sample_data(pseq_tmp2_cox)$DEATH <- death_vars$DEATH
  sample_data(pseq_tmp2_cox)$DEATH_AGEDIFF <- death_vars$DEATH_AGEDIFF

  cox_data <- t(abundances(pseq_tmp2_cox))

  cox_data <- cbind(cox_data,
                          meta(pseq_tmp2_cox)[, c(covariates, "DEATH", "DEATH_AGEDIFF")])
  cox_data
}