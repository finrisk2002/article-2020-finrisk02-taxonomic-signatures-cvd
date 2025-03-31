# cause_of_death_enterobacteriaceae.R
Cox regression and forest plot for the hazard ratio of Ensterobacteriaceae group in various causes of death

# create_figure_genus_panel_srf_importance.R
Extended Data Fig. 5: Examples of hazard ratios for genera and random forest importance bar panel

# create_figure_PC3_drivers_east_west.R
Extended Data Fig. 4: PC3 drivers in east and west

# create_figure_PC3_east_west_mortality.R
Extended Data Fig. 2: HR for death vs. PC3 in east and west

# create_figure_PC_drivers.R
Extended Data Fig. 3: PC 1-3 drivers

# create_figure_pca_mortality_association.R
Figure 2.: HR for death vs. PC 1-3

# functional_cox.R
Cox regression: functional metabolites

# genus_cox.R
Cox regression: individual genera. 

# misc.checks.R
Minor miscellaneous checks, mostly for review letter

# PC3_supplementary_checks.R
Mostly additional checks concerning PC3 hazard for the review letter:
Cox mortality ~ covariates + PC3 with..
1) Prints results for all covariates
2) Including Healthy Diet Index as covariate
3) Excluding prevalent cancers
4) Hazard ratio between Q1 and Q4
5) 

# PCA_cox.R
Cox regression with predictors: PC axes 1-3, Shannon diversity and Observed rihcness. 

# subnet_cox.R
Cox regression: subnet abundances

# survival_main.R
Running this script runs the entire survival analysis, creates the figures and table.

# survival_random_forest.R
Executes survival random forest with predictors: core genera, covariates, core and covariates

# survival_random_forest_cross_validation.R
5-fold cross validation on survival random forest. Predictor sets used: core genera, covariates, core and covariates. Also computes Harrell's c-statistics and compares these (paired t-test) between the predictor sets. 