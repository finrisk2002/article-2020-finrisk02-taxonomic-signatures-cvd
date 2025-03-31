# Common R functions


## Cox ********************************************* ####

## I'm not quite sure if all of these functions are used, but better to keep all here for now (VL)

## Plan would be to move major functions to their own files at R_functions/

## Generate formula, fit Cox PH and wrangle results. 
# data = data.frame with all needed variables as columns
# predictors = vector of variables whose effect we are interested in
# covariates = vector of controlled covariates
# status = 0/1 corresponding to no event/event; alive/dead etc.
# time_to_event = follow-up time until event or censoring
# splines = TRUE/FALSE; should the effect of predictor be modelled non-linearly? 
# If TRUE, both linear and non-linear models are fit and for each predictor the fit with lower p-value is chosen.
# normalize = TRUE/FALSE; standardize (x-mean(x)/sd(x)) predictor before fitting model?
# This makes comparison of hazard ratios easier
# test_ph_assumption = TRUE/FALSE; if TRUE the function also returns a data frame with variable pairs that potentially
# do not meet the assumption of proportional hazards.


cox_regression <- function(data,
                           predictors,
                           covariates,
                           status = "DEATH",
                           time_to_event = "DEATH_AGEDIFF",
                           splines,
                           normalize) {
  
  
  if(normalize) {
    
    if(class(data[, predictors]) == "numeric") {
      
      x <- data[, predictors]
      data[, predictors] <- (x - mean(x, na.rm = T))/sd(x, na.rm = T) 
    } else {
      data[, predictors] <- apply(data[, predictors], 2, FUN = function(x) {(x - mean(x, na.rm = T))/sd(x, na.rm = T) })
    }
    
  }
  
  ## Formulas ***************************
  linear_formulas <- lapply(predictors, function(x) {
    formula_data <- deparse(substitute(data))
    
    formula <- paste0("Surv(", formula_data, "$", time_to_event, ", ",formula_data,"$", status, ") ~ ",paste(covariates, collapse = "+"), " + ",x)
    
    return(formula)
    
  }) %>% 
    set_names(predictors)
  
  if(splines) {
    spline_formulas <- lapply(predictors, function(x) {
      formula_data <- deparse(substitute(data))
      
      formula <- paste0("Surv(", formula_data, "$", time_to_event, ", ",formula_data,"$", status, ") ~ ",paste(covariates, collapse = "+"), " + pspline(",x, ")")
      
      return(formula)
      
    }) %>% 
      set_names(predictors)
  }
  
  ## Cox regression *********************
  print("Cox")
  linear_cox_fit <- lapply(linear_formulas, function(x) {
    
    coxph(as.formula(x), data=data, x=TRUE)
    
  })
  
  if(splines) {
    spline_cox_fit <- lapply(spline_formulas, function(x) {
      
      coxph(as.formula(x), data=data, x=TRUE)
      
    })
    
    return(list(linear_cox_fit = linear_cox_fit, spline_cox_fit = spline_cox_fit))
  } else {
    
    
    return(linear_cox_fit = linear_cox_fit)
  }
  
  
  
}


cox_results <- function(data,
                        predictors,
                        covariates,
                        status = "DEATH",
                        time_to_event = "DEATH_AGEDIFF",
                        alpha_level,
                        splines, 
                        fit_list) {
  
  
  linear_cox_fit <- fit_list[["linear_cox_fit"]]
  if(splines) {spline_cox_fit <- fit_list[["spline_cox_fit"]]}
  
  
  ## Results *****************************
  print("Results")
  results <- lapply(predictors, function(x) {
    
    df <- summary(linear_cox_fit[[x]])$coefficients %>% as.data.frame()
    df <- df[nrow(df), ] %>% 
      select(coef, "se(coef)", "z", "Pr(>|z|)") %>% 
      set_colnames(c("coef", "se_coef", "test_statistic_value", "p")) %>% 
      mutate(test_statistic = "Wald")
    
    
    
    if(splines) {
      spline_df <- summary(spline_cox_fit[[x]])$coefficients %>%
        as.data.frame() %>% 
        select(coef, "se(coef)", "Chisq", p) %>% 
        set_colnames(c("coef", "se_coef", "test_statistic_value", "p")) %>% 
        mutate(test_statistic = "Chisq") %>% 
        slice(nrow(.))
      
      # filter less significant association
      df <- rbind(df, spline_df) %>% 
        mutate(association = c("linear", "non-linear")) 
    } 
    
    
    
    df <- df %>% 
      mutate(predictor = x)
    
    
    
    
  }) %>% 
    do.call(rbind, .)
  
  # Multiple testing correction
  results <- results %>%
    group_by(association) %>%
    mutate(P_adjusted = p.adjust(p, "BH")) %>%
    ungroup() %>% 
    group_by(predictor) 
  
  
  # Results in neat form for presentation
  neat_results <- results %>% 
    # filter(p == min(p)) %>% 
    ungroup() %>% 
    mutate(HR = round(exp(coef),3)) %>% 
    mutate(HR_lower_95 = round(exp(coef - 1.96*se_coef), 3),
           HR_upper_95 = round(exp(coef + 1.96*se_coef), 3),
           P = round(p, 5),
           Coefficient = round(coef, 3),
           "Coefficient SE" = round(se_coef, 3)) %>% 
    mutate(HR = paste0(HR, " (95% CI, ", HR_lower_95, "-", HR_upper_95, ")")) %>% 
    select(Predictor =predictor, Coefficient, "Coefficient SE", HR, "P_adjusted", test_statistic_value, test_statistic) %>%
    mutate(HR = ifelse(is.na(Coefficient), NA, HR), Association = ifelse(is.na(Coefficient), "non-linear", "linear"))  %>% 
    filter(P_adjusted < alpha_level) %>% 
    arrange(Association, P_adjusted) %>% 
    set_colnames(c("Predictor", "Coefficient", "Coefficient SE", "HR", "P (adjusted)", "Test Statistic Value", "Test Statistic", "Association"))
  
  
  # Results in a form more convenient for further manipulations
  results <- results %>% 
    ungroup %>% 
    mutate(PH = exp(coef)) %>% 
    mutate(p_adj = P_adjusted) %>% 
    mutate(linearity = ifelse(is.na(coef), "non-linear", "linear"),
           posneg = ifelse(is.na(coef), "nonlinear",
                           ifelse(coef > 0, "bad", "good")))
  
  linear_results <- results %>% filter(linearity == "linear")
  non_linear_results <- results %>% filter(linearity == "non-linear")
  
  if(nrow(neat_results) == 0) {
    return(list(results = results,
                linear_results = linear_results,
                non_linear_results = non_linear_results))
    
  }
  
  return(list(neat_results = neat_results, 
              results = results,
              linear_results = linear_results,
              non_linear_results = non_linear_results))
  
  
}

cox_ph_assumptions <- function(data,
                               predictors,
                               covariates,
                               status = "DEATH",
                               time_to_event = "DEATH_AGEDIFF",
                               splines,
                               fit_list) {
  
  
  linear_cox_fit <- fit_list[["linear_cox_fit"]]
  if(splines) {spline_cox_fit <- fit_list[["spline_cox_fit"]]}
  
  
  ph_assumption <-  lapply(predictors, function(m) {
    test <- cox.zph(linear_cox_fit[[m]])
    
    p_values <- test$table[, "p"]
    
    # significant cases
    x <- which(p_values < 1)
    
    if(length(x) == 0) {
      return(NULL)
    }
    
    df <- data.frame(feature = m, variable_not_ph = names(x), p_value = p_values[x])
    
  }) %>% 
    do.call(rbind, .) %>%
    mutate(p_adj = p.adjust(p_value, "BH")) %>% 
    filter(p_value < alpha_level)
  
  
  return(ph_assumption)
  
}


# Wrap the above three functions into one
cox_wrapper <- function(data,
                        predictors,
                        covariates,
                        status = "DEATH",
                        time_to_event = "DEATH_AGEDIFF",
                        alpha_level = alpha_level,
                        splines,
                        normalize,
                        test_ph_assumption = FALSE) {
  
  
  if(normalize) {
    
    if(class(data[, predictors]) == "numeric") {
      
      x <- data[, predictors]
      data[, predictors] <- (x - mean(x, na.rm = T))/sd(x, na.rm = T) 
    } else {
      data[, predictors] <- apply(data[, predictors], 2, FUN = function(x) {(x - mean(x, na.rm = T))/sd(x, na.rm = T) })
    }
    
  }
  
  ## Formulas ***************************
  linear_formulas <- lapply(predictors, function(x) {
    formula_data <- deparse(substitute(data))
    
    formula <- paste0("Surv(", formula_data, "$", time_to_event, ", ",formula_data,"$", status, ") ~ ",paste(covariates, collapse = "+"), " + ",x)
    
    return(formula)
    
  }) %>% 
    set_names(predictors)
  
  if(splines) {
    spline_formulas <- lapply(predictors, function(x) {
      formula_data <- deparse(substitute(data))
      
      formula <- paste0("Surv(", formula_data, "$", time_to_event, ", ",formula_data,"$", status, ") ~ ",paste(covariates, collapse = "+"), " + pspline(",x, ")")
      
      return(formula)
      
    }) %>% 
      set_names(predictors)
  }
  
  ## Cox regression *********************
  print("Cox")
  linear_cox_fit <- lapply(linear_formulas, function(x) {
    
    coxph(as.formula(x), data=data, x=TRUE)
    
  })
  
  if(splines) {
    spline_cox_fit <- lapply(spline_formulas, function(x) {
      
      coxph(as.formula(x), data=data, x=TRUE)
      
    })
  }
  
  
  ## Check PH assumptions ****************
  if(test_ph_assumption) {
    print("PH assumptions")
    ph_assumption <-  lapply(predictors, function(m) {
      test <- cox.zph(linear_cox_fit[[m]])
      
      p_values <- test$table[, "p"]
      
      # significant cases
      x <- which(p_values < 1)
      
      if(length(x) == 0) {
        return(NULL)
      }
      
      df <- data.frame(feature = m, variable_not_ph = names(x), p_value = p_values[x])
      
    }) %>% 
      do.call(rbind, .) %>%
      mutate(p_adj = p.adjust(p_value, "BH")) %>% 
      filter(p_value < alpha_level)
    
  }
  
  ## Results *****************************
  print("Results")
  results <- lapply(predictors, function(x) {
    
    df <- summary(linear_cox_fit[[x]])$coefficients %>% as.data.frame()
    df <- df[nrow(df), ] %>% 
      select(coef, "se(coef)", "z", "Pr(>|z|)") %>% 
      set_colnames(c("coef", "se_coef", "test_stat_value", "p")) %>% 
      mutate(test_stat = "Wald")
    
    
    
    if(splines) {
      spline_df <- summary(spline_cox_fit[[x]])$coefficients %>%
        as.data.frame() %>% 
        select(coef, "se(coef)", Chisq,  p) %>% 
        set_colnames(c("coef", "se_coef", "test_stat_value", "p")) %>% 
        mutate(test_stat = "Chisq") %>% 
        slice(nrow(.))
      
      df <- rbind(df, spline_df) %>% 
        mutate(association = c("linear", "non-linear")) 
    } 
    
    
    
    df <- df %>% 
      mutate(predictor = x)
    
    
    
    
  }) %>% 
    do.call(rbind, .) 
  
  # Multiple testing correction
  results <- results %>%
    group_by(association) %>% 
    mutate(P_adjusted = p.adjust(p, "BH")) %>% 
    ungroup() %>% 
    group_by(predictor) 
  
  
  
  
  
  
  
  
  # Results in neat form for presentation
  neat_results <- results %>% 
    # filter(p == min(p)) %>% 
    ungroup() %>% 
    mutate(HR = round(exp(coef),3)) %>% 
    mutate(HR_lower_95 = round(exp(coef - 1.96*se_coef), 3),
           HR_upper_95 = round(exp(coef + 1.96*se_coef), 3),
           P = round(p, 5),
           Coefficient = round(coef, 3),
           "Coefficient SE" = round(se_coef, 3)) %>% 
    mutate(HR = paste0(HR, " (95% CI, ", HR_lower_95, "-", HR_upper_95, ")")) %>% 
    select(Predictor = predictor, Coefficient, "Coefficient SE", HR, "P_adjusted", "test_stat_value", "test_stat") %>%
    mutate(HR = ifelse(is.na(Coefficient), NA, HR), Association = ifelse(is.na(Coefficient), "non-linear", "linear"))  %>% 
    filter(P_adjusted < alpha_level) %>% 
    arrange(Association, P_adjusted) %>% 
    set_colnames(c("Predictor", "Coefficient", "Coefficient SE", "HR", "P (adjusted)", "Test Statistic Value", "Test Statistic", "Association"))
  
  
  # Results in a form more convenient for further manipulations
  results <- results %>% 
    ungroup %>% 
    mutate(PH = exp(coef)) %>% 
    mutate(p_adj = P_adjusted) %>% 
    mutate(linearity = ifelse(is.na(coef), "non-linear", "linear"),
           posneg = ifelse(is.na(coef), "nonlinear",
                           ifelse(coef < 0, "bad", "good")))
  
  
  
  
  
  
  if(nrow(neat_results) == 0) {
    return(list(results = results))
  }
  
  if(test_ph_assumption) {
    
    if(nrow(neat_results) == 0) {
      return(list(results = results,
                  ph_assumption = ph_assumption))
    }
    
    return(list(neat_results = neat_results, 
                results = results,
                ph_assumption = ph_assumption))
  }
  
  return(list(neat_results = neat_results, 
              results = results))
  
}


# wrapper2: returns covariate results as well 
cox_wrapper_with_covariates <- function(data,
                                        predictors ,
                                        covariates,
                                        status = "DEATH",
                                        time_to_event = "DEATH_AGEDIFF",
                                        alpha_level = alpha_level,
                                        splines,
                                        normalize, 
                                        p_adjust = TRUE) {
  
  
  if(normalize) {
    
    if(class(data[, predictors]) == "numeric") {
      
      x <- data[, predictors]
      data[, predictors] <- (x - mean(x, na.rm = T))/sd(x, na.rm = T) 
    } else {
      data[, predictors] <- apply(data[, predictors],
                                  2,
                                  FUN = function(x) {
                                    (x - mean(x, na.rm = T))/sd(x, na.rm = T)
                                  })
    }
    
    
    for(co in covariates) {
      
      if(class(data[, co]) == "numeric") {
        
        x <- data[, co]
        data[, co] <- (x - mean(x, na.rm = T))/sd(x, na.rm = T) 
      }
      
    }
    
  }
  
  ## Formulas ***************************
  linear_formulas <- lapply(predictors, function(x) {
    formula_data <- deparse(substitute(data))
    
    formula <- paste0("Surv(", formula_data, "$", time_to_event, ", ",formula_data,"$", status, ") ~ ",paste(covariates, collapse = "+"), " + ",x)
    
    return(formula)
    
  }) %>% 
    set_names(predictors)
  
  if(splines) {
    spline_formulas <- lapply(predictors, function(x) {
      formula_data <- deparse(substitute(data))
      
      formula <- paste0("Surv(", formula_data, "$", time_to_event, ", ",formula_data,"$", status, ") ~ ",paste(covariates, collapse = "+"), " + pspline(",x, ")")
      
      return(formula)
      
    }) %>% 
      set_names(predictors)
  }
  
  ## Cox regression *********************
  print("Cox")
  linear_cox_fit <- lapply(linear_formulas, function(x) {
    
    coxph(as.formula(x), data=data, x=TRUE)
    
  })
  
  if(splines) {
    spline_cox_fit <- lapply(spline_formulas, function(x) {
      
      coxph(as.formula(x), data=data, x=TRUE)
      
    })
  }
  
  
  ## Results *****************************
  print("Results")
  results <- lapply(predictors, function(x) {
    
    df <- summary(linear_cox_fit[[x]])$coefficients %>% as.data.frame()
    df <- df %>% 
      rownames_to_column(var = "variable") %>% 
      select(variable, coef, "se(coef)", "z", "Pr(>|z|)") %>%
      set_colnames(c("variable", "coef", "se_coef", "test_stat_value", "p")) %>%
      mutate(test_stat = "Wald")
    
    
    
    # if(splines) {
    #   spline_df <- summary(spline_cox_fit[[x]])$coefficients %>%
    #     as.data.frame() %>% 
    #     rownames_to_column(var = "variable") %>% 
    #     select(variable, coef, "se(coef)", Chisq,  p) %>% 
    #     set_colnames(c("variable", "coef", "se_coef", "test_stat_value", "p")) %>% 
    #     mutate(test_stat = "Chisq") %>% 
    #     slice(-(nrow(.) - 1))
    #   
    #   df <- rbind(df, spline_df) %>% 
    #     mutate(association = c("linear", "non-linear")) 
    # } 
    
    
    
    return(df)
    
  }) %>% 
    do.call(rbind, .) 
  
  # Multiple testing correction
  # results <- results %>%
  # mutate(P_adjusted = p.adjust(p, "BH"))
  
  
  
  
  
  
  
  
  
  
  
  # Results in neat form for presentation
  neat_results <- results %>% 
    # filter(p == min(p)) %>% 
    ungroup() %>% 
    mutate(HR = round(exp(coef),3)) %>% 
    mutate(HR_lower_95 = round(exp(coef - 1.96*se_coef), 3),
           HR_upper_95 = round(exp(coef + 1.96*se_coef), 3),
           # P = round(p, 5),
           P = p,
           Coefficient = round(coef, 3),
           "Coefficient SE" = round(se_coef, 3), 
           Variable = variable) %>% 
    mutate(HR = paste0(HR, " (95% CI, ", HR_lower_95, "-", HR_upper_95, ")")) %>% 
    select(Variable, Coefficient, "Coefficient SE", HR, "P", "test_stat_value", "test_stat") %>%
    mutate(HR = ifelse(is.na(Coefficient), NA, HR), Association = ifelse(is.na(Coefficient), "non-linear", "linear"))  %>% 
    # filter(P_adjusted < alpha_level) %>%
    arrange(Association, P) %>% 
    set_colnames(c("Variable", "Coefficient", "Coefficient SE", "HR", "P", "Test Statistic Value", "Test Statistic", "Association"))
  
  
  # Results in a form more convenient for further manipulations
  results <- results %>% 
    ungroup %>% 
    mutate(PH = exp(coef)) %>% 
    mutate(linearity = ifelse(is.na(coef), "non-linear", "linear"))
  
  
  
  
  
  
  if(nrow(neat_results) == 0) {
    return(list(results = results))
  }
  
  if(test_ph_assumption) {
    
    if(nrow(neat_results) == 0) {
      return(list(results = results,
                  ph_assumption = ph_assumption))
    }
    
    return(list(neat_results = neat_results, 
                results = results,
                ph_assumption = ph_assumption))
  }
  
  return(list(neat_results = neat_results, 
              results = results))
  
}

# Wrapper for categorical predictor. Compares last level to the first
cox_wrapper_categorical <- function(data,
                                    predictors ,
                                    covariates,
                                    status = "DEATH",
                                    time_to_event = "DEATH_AGEDIFF",
                                    alpha_level = alpha_level,
                                    normalize) {
  
  
  if(normalize) {
    
    if(class(data[, predictors]) == "numeric") {
      
      x <- data[, predictors]
      data[, predictors] <- (x - mean(x, na.rm = T))/sd(x, na.rm = T) 
    }
    
    for(co in covariates) {
      
      if(class(data[, co]) == "numeric") {
        
        x <- data[, co]
        data[, co] <- (x - mean(x, na.rm = T))/sd(x, na.rm = T) 
      }
      
    }
    
  }
  
  ## Formulas ***************************
  print("formulas")
  linear_formulas <- lapply(predictors, function(x) {
    formula_data <- deparse(substitute(data))
    
    formula <- paste0("Surv(", formula_data, "$", time_to_event, ", ",formula_data,"$", status, ") ~ ",paste(covariates, collapse = "+"), " + ",x)
    
    return(formula)
    
  }) %>% 
    set_names(predictors)
  
  
  ## Cox regression *********************
  print("Cox")
  linear_cox_fit <- lapply(linear_formulas, function(x) {
    
    coxph(as.formula(x), data=data, x=TRUE)
    
  })
  
  
  ## Results *****************************
  print("Results")
  results <- lapply(predictors, function(x) {
    
    df <- summary(linear_cox_fit[[x]])$coefficients %>% as.data.frame()
    df <- df %>% 
      rownames_to_column(var = "variable") %>% 
      select(variable, coef, "se(coef)", "z", "Pr(>|z|)") %>%
      set_colnames(c("variable", "coef", "se_coef", "test_stat_value", "p")) %>%
      mutate(test_stat = "Wald")
    
    
    
    # if(splines) {
    #   spline_df <- summary(spline_cox_fit[[x]])$coefficients %>%
    #     as.data.frame() %>% 
    #     rownames_to_column(var = "variable") %>% 
    #     select(variable, coef, "se(coef)", Chisq,  p) %>% 
    #     set_colnames(c("variable", "coef", "se_coef", "test_stat_value", "p")) %>% 
    #     mutate(test_stat = "Chisq") %>% 
    #     slice(-(nrow(.) - 1))
    #   
    #   df <- rbind(df, spline_df) %>% 
    #     mutate(association = c("linear", "non-linear")) 
    # } 
    
    
    
    return(df)
    
  }) %>% 
    do.call(rbind, .) 
  
  # Results in neat form for presentation
  print("neat_resiults")
  neat_results <- results %>% 
    # filter(p == min(p)) %>% 
    ungroup() %>% 
    mutate(HR = round(exp(coef),3)) %>% 
    mutate(HR_lower_95 = round(exp(coef - 1.96*se_coef), 3),
           HR_upper_95 = round(exp(coef + 1.96*se_coef), 3),
           # P = round(p, 5),
           P = p,
           Coefficient = round(coef, 3),
           "Coefficient SE" = round(se_coef, 3), 
           Variable = variable) %>% 
    mutate(HR = paste0(HR, " (95% CI, ", HR_lower_95, "-", HR_upper_95, ")")) %>% 
    select(Variable, Coefficient, "Coefficient SE", HR, "P", "test_stat_value", "test_stat") %>%
    mutate(HR = ifelse(is.na(Coefficient), NA, HR), Association = ifelse(is.na(Coefficient), "non-linear", "linear"))  %>% 
    # filter(P_adjusted < alpha_level) %>%
    arrange(Association, P) %>% 
    set_colnames(c("Variable", "Coefficient", "Coefficient SE", "HR", "P", "Test Statistic Value", "Test Statistic", "Association"))
  
  
  # Results in a form more convenient for further manipulations
  print("results")
  results <- results %>% 
    ungroup %>% 
    mutate(PH = exp(coef)) %>% 
    mutate(linearity = ifelse(is.na(coef), "non-linear", "linear"))
  
  
  
  
  print("retursn")
  
  if(nrow(neat_results) == 0) {
    return(list(results = results))
  }
  
  
  return(list(neat_results = neat_results, 
              results = results))
  
}


# Hazard ratio plotter
# Note that the fit needs to be explicitely written and not from a loop.
plot_HR <- function(fit,
                    pred,
                    data=microbiome::meta(pseq),
                    ref_value=NULL,
                    rects = FALSE,
                    ribbon_color = "#5C88DA",
                    title = TRUE, 
                    return_data = FALSE) {
  
  
  # clean predicator name
  pred_name <- clean_genus_names(pred)
  
  # rename data
  mydata <- data
  
  # not sure if this is redundant
  hr <- smoothHR(data=(mydata)[, c(covariates, pred, "DEATH_AGEDIFF", "DEATH")], coxfit = fit)
  
  # use predictor median as reference value, unless defined in function call
  if(is.null(ref_value)) {
    ref_value <- median(mydata[, pred])
  }
  
  # range of predictor for graph
  r <- range(hr$dataset[, pred])
  
  
  
  # simulate values and edit data frame
  fit_sim <- predict(hr,
                     predictor = pred,
                     pred.value = ref_value,
                     prediction.values = seq(from=r[1], to=r[2], length.out=100), data=mydata)
  
  fit_sim_df <- fit_sim %>% as.data.frame()
  colnames(fit_sim_df) <- c("pred", "ln_hr", "ln_lower95", "ln_upper95")
  
  
  
  # natural logarithm
  fit_sim_df <- fit_sim_df %>%
    mutate(hr = exp(ln_hr), lower95=exp(ln_lower95), upper95=exp(ln_upper95)) %>%
    mutate(log2_hr=log2(hr), log2_lower95=log2(lower95), log2_upper95=log2(upper95))
  
  if(return_data) {
    return(fit_sim_df)
  }
  
  # plot
  if(!isTRUE(title)) {
    p <- ggplot(data.frame(fit_sim_df), aes(x=pred, y=log2_hr)) +
      geom_line(size = 1) +
      geom_hline(yintercept = 0, linetype="dashed") +
      geom_ribbon(aes(ymin = fit_sim_df$log2_lower95, ymax = fit_sim_df$log2_upper95),
                  linetype = 2,
                  alpha = .5,
                  fill = ribbon_color) +
      scale_y_continuous(breaks=-2:6, labels = 2^(-2:6)) +
      labs(y="HR for Death") + 
      coord_cartesian(xlim = c(r[1], r[2]))
    
    
    
  } else {
    
    p <- ggplot(data.frame(fit_sim_df), aes(x=pred, y=log2_hr)) +
      geom_line(size = 1) +
      geom_hline(yintercept = 0, linetype="dashed") +
      geom_ribbon(aes(ymin = fit_sim_df$log2_lower95, ymax = fit_sim_df$log2_upper95),
                  linetype = 2,
                  alpha = .5,
                  fill = ribbon_color) +
      scale_y_continuous(breaks=-2:6, labels = 2^(-2:6)) +
      labs(y="HR for Death", title=pred_name) + 
      coord_cartesian(xlim = c(r[1], r[2]))
  }
  
  
  
  # Add quantile regions  
  if(rects) {
    df <- mydata %>% select(one_of(pred))
    
    q <- quantile(df[, pred], probs = c(0, .2, .4, .6, .8, 1))
    
    rects <- data.frame(xstart = c(q[1], q[2], q[3], q[4], q[5]), xend = c(q[2], q[3], q[4], q[5], q[6]), quantile = as.character(1:5))
    
    p2 <-  ggplot() + geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = quantile), alpha = 0.4) + scale_fill_manual(values=c('#fef0d9','#fdcc8a','#fc8d59','#e34a33','#b30000')) +
      geom_line(data = fit_sim_df, aes(x=pred, y=log2_hr)) +
      geom_hline(yintercept = 0, linetype="dashed") +
      geom_ribbon(data = fit_sim_df,aes(x=pred, ymin=fit_sim_df$log2_lower95, ymax=fit_sim_df$log2_upper95), linetype=2, alpha=0.1) +
      scale_y_continuous(breaks=-2:6, labels = 2^(-2:6), limits = c(-3,4)) +
      # scale_x_continuous(breaks=-5:12, limits = c(r[1], r[2])) +
      labs(y="Proportional hazard", x=pred_name)+ theme_bw()  
    
    return(p2)
  }
  
  return(p)
  
}







# Generate Cox model formulas
cox_formula <- function(data = meta, predictor, covariates, splines = TRUE) {
  
  pred <- ifelse(class(data[, predictor])=="numeric" & splines == TRUE, paste0("+ pspline(", predictor,")"), paste("+",predictor)) 
  
  formula_data <- deparse(substitute(data))
  
  formula <- paste0("Surv(",formula_data,"$DEATH_AGEDIFF, ",formula_data,"$Death) ~",paste(covariates, collapse = "+"),pred)
  
  return(formula)
  
}


functional_cox_formula <- function(data, predictor, covariates) {
  
  pred <- ifelse(class(data[, predictor])=="numeric", paste0("+ pspline(", predictor,")"), paste("+",predictor)) 
  
  formula <- paste("Surv(meta$DEATH_AGEDIFF, meta$Death) ~",paste(covariates, collapse = "+"),pred)
  
  return(formula)
  
}


# Get results
get_cox_results <- function(var, model_list, data=NULL) {
  
  if(!is.null(data)) {
    meta <- data
  }
  
  # separate cases for numeric and factor variables
  if(class(meta[, var])=="numeric") {
    
    # number of rows in the summary. 
    # Necessary to extract rows form summary by index as
    # long predictor names can be cut
    
    
    df <- summary(model_list[[var]])$coefficients %>% 
      as.data.frame %>% 
      rownames_to_column(var = "predictor")
    
    df %>% 
      filter(predictor == var) %>% 
      select(predictor, coef, "se(coef)", "Pr(>|z|)") %>%
      set_colnames(c("predictor", "coef", "se_coef", "p"))
    
  } else if(class(meta[, var])=="factor") {
    
    summary(model_list[[var]])$coefficients %>% 
      as.data.frame %>% 
      rownames_to_column(var = "predictor") %>%
      filter(grepl(var, predictor)) %>%
      select("Pr(>|z|)", coef, "se(coef)") %>%
      set_colnames(c("p", "coef", "se_coef")) %>% 
      mutate(predictor=var, linearity=NA, predictor=var, pred_class="factor")
    
  } else warning("Check variable class")
  
  
  
}



# Cox regression cross validatation. Computes the ROC and AUC. 
cox_cv <- function(data,
                   predictor,
                   response_var, 
                   response_time = "DEATH_AGEDIFF",
                   covariates = covariates,
                   n_folds = 5,
                   seed = 11235) {
  
  # set seed for reproducibility
  set.seed(seed)
  
  # Data **************************************************************************
  
  # mydata is data with used variables filtered
  mydata <- data[, c(response_var,
                     response_time,
                     predictor,
                     covariates)] %>% 
    drop_na()
  
  
  # modify variable names so upcoming function calls will work
  mydata$response_time = mydata[, response_time]
  mydata[, response_time] <- NULL
  mydata$response_var = mydata[, response_var]
  mydata[, response_var] <- NULL
  
  
  
  # Folds *************************************************************************
  # Randomly shuffle the data
  mydata <- mydata[sample(nrow(mydata)),]
  
  # Create n_folds equally sized folds
  folds <- cut(seq(1, nrow(mydata)), breaks = n_folds, labels = FALSE)
  
  
  # Cross validate ***************************************************************
  print("CV")
  cox_cv <- lapply(1:n_folds, function(fold) {
    print(fold)
    
    # training samples
    train <- folds != fold
    
    # fit cox
    
    formula <- paste0("Surv(response_time, response_var) ~ BL_AGE + BMI + MEN + pspline(",predictor,")") %>% 
      as.formula()
    
    cox_train <- coxph(formula = formula, 
                       data = mydata[train, ],
                       x = TRUE)
    
    # Predict in test set
    cox_pred <- predict(cox_train, newdata = mydata[!train, ])
    
    return(cox_pred)
    
  }) 
  
  
  # ROC & AUC ********************************************************************
  print("ROC & AUC")
  ROC_df <- lapply(1:n_folds, function(fold) {
    print(fold)
    
    train <- folds != fold
    
    # Combine
    dat <- cbind(mydata[!train, ], survival_pred = cox_cv[[fold]])
    
    roc_res <- with(dat,
                    survivalROC(Stime        = response_time,
                                status       = response_var,
                                marker       = survival_pred,
                                predict.time = max(mydata[, "response_time"]),
                                method       = "KM"))       # KM method without smoothing
    
    df <- data.frame(predictor = predictor,
                     FPR = roc_res[["FP"]],
                     TPR = roc_res[["TP"]],
                     fold = rep(fold, length(roc_res[["TP"]])),
                     AUC = rep(roc_res[["AUC"]],length(roc_res[["TP"]])))
    return(df)
    
  }) %>% do.call(rbind, .)
  
  return(ROC_df)
}

## Survival Random Forest ************************** ####

# Compute the harrell's c in n_folds cross-validation.
random_survival_forest_cv <- function(data,
                                       predictor_list,
                                       response_var,
                                       response_time = "DEATH_AGEDIFF",
                                       n_folds = 5,
                                       seed = 11235,
                                       save_file = "") {
  
  
  
  
  # set seed for reproducibility
  set.seed(seed)
  
  # Data **************************************************************************
  # save results to 
  survival_RF <- list()
  
  # mydata is data with used variables filtered
  mydata <- data[, c(response_var,
                     response_time,
                     predictor_list %>% unlist %>% unique())] %>% 
    drop_na()
  
  # drop subjects with no response time
  mydata <- mydata %>% drop_na(response_time)
  
  # modify variable names so upcoming function calls will work
  mydata$response_time = mydata[, response_time]
  mydata[, response_time] <- NULL
  mydata$response_var = mydata[, response_var]
  mydata[, response_var] <- NULL
  
  
  # Folds *************************************************************************
  # Randomly shuffle the data
  mydata <- mydata[sample(nrow(mydata)),]
  
  # Create n_folds equally sized folds
  folds <- cut(seq(1, nrow(mydata)), breaks = n_folds, labels = FALSE)
  
  
  # Loop over list of predictors. Use for() to save intermediate steps and ease debugging
  for(subs in names(predictor_list)) { 
    print(subs)
    
    # subset data 
    rf_data <- mydata[, c("response_time", "response_var", predictor_list[[subs]])]
    
    
    # Cross validate ***************************************************************
    print("CV")
    rf_cv <- lapply(1:n_folds, function(fold) {
      print(fold)
      
      # training samples
      train <- folds != fold
      
      # grow random forest
      forest <- rfsrc(Surv(response_time, response_var) ~ .,
                      data = rf_data[train, ],
                      importance = TRUE)
      
      # predict in test set
      pred <- predict(forest,
                      newdata = rf_data[!train, ])
      
      return(list(forest_importance = forest$importance,
                  prediction_values = pred$predicted,
                  train_samples = train,
                  harrell_c = 1 - forest$err.rate[length(forest$err.rate)]))
      
    })
    
    # ROC & AUC ********************************************************************
    print("ROC & AUC")
    ROC_df <- lapply(1:n_folds, function(fold) {
      print(fold)
      train <- folds != fold
      
      roc_data <- cbind(rf_data[!train, ], survival_pred = rf_cv[[fold]]$prediction_values)
      roc_res <- with(roc_data,
                      survivalROC(Stime = response_time,
                                  status = response_var,
                                  marker = survival_pred,
                                  predict.time = max(mydata[, "response_time"]),
                                  method = "KM"))       # KM method without smoothing
      
      df <- data.frame(predictor = subs,
                       FPR = roc_res[["FP"]],
                       TPR = roc_res[["TP"]],
                       fold = rep(fold, length(roc_res[["TP"]])),
                       AUC = rep(roc_res[["AUC"]],length(roc_res[["TP"]])))
      return(df)
      
    }) %>% do.call(rbind, .)
    
    
    # Results to list
    survival_RF[[subs]] <- list(cv = rf_cv, roc = ROC_df)
    
  }
  
  return(survival_RF)
  
}







## MISC ******************************************** ####


get_causes_of_death <- function(x) {
  
  if(class(x) == "phyloseq") {
    if (!("K_TPKS" %in% sample_variables(x)) && "UCOD" %in% sample_variables(x)) {
      warning("no raw causes of death variable K_TPKS in phyloseq, but UCOD present. using get_causes_of_death_ucod instead")
      return(get_causes_of_death_ucod(x))
    }
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

get_causes_of_death_ucod <- function(x) {
  
  if(class(x) == "phyloseq") {
    if ("UCOD" %in% sample_variables(x)) {
      message("UCOD present, renaming it cause_of_death")
      cause_of_death <- meta(x)[, "UCOD"]
      sample_data(x) <- cbind(meta(x), cause_of_death)
      return(x)
    } else {
      stop("UCOD not found in phyloseq sample_data. Are you looking for get_causes_of_death()?")
    }
  }
  stop("not a phyloseq. Are you looking for get_causes_of_death()?")
}


# Phyloseq editor 
# Takes a phyloseq object and adds pca components, diversity measures, edit taxa names

add_many_things_to_meta_data <- function(pseq, species_pseq) {
  
  my_pseq <- pseq
  
  
  ## Diversities
  my_pseq <- add_diversities_to_meta_data(my_pseq)
  
  
  ## Taxa names: Bacteroides (Bacteria) --> Bacteroides_Bacteria
  taxa_names(my_pseq) <- dirty_genus_names(taxa_names(my_pseq))
  
  
  ## CLR transformed abundances
  my_pseq <- add_clr_abundances_to_meta_data(my_pseq)
  
  
  ## Scaled PCA components, continuous
  my_pseq <- add_pca_components_to_meta_data(my_pseq, species_pseq)
  
  
  return(my_pseq)
  
}

add_pca_components_to_meta_data <- function(pseq, species_pseq, n_axes = 3, check_existence = TRUE) {
  
  set.seed(12345)
  
  my_pseq <- pseq
  
  # Check if PCA output file exists. If not, then compute
  if(file.exists("output/raw_output/species_pca.rds") & check_existence) {
    pca_clr <- readRDS(file = "output/raw_output/species_pca.rds")
  } else {
    
    print("Doing PCA. This may take a while!")
    
    pca_clr <- prcomp(microbiome::transform(
      microbiome::transform(species_pseq,
                            "compositional"),
      "clr"
    ) %>%
      abundances %>%
      t)
    
    # Flip PCs to get a consistent directions
    # (sign of PCA axis weight is arbitrary:
    # but we noticed that directions might flip randomly even with constant seed)
    # For consistency, Candidatus Korarchaeota should have negative PC1, positive PC2, positive PC3
    if (pca_clr$rotation[1, "PC1"] > 0) {
      pca_clr$rotation[, "PC1"] <- -1*pca_clr$rotation[, "PC1"]
      pca_clr$x[, "PC1"] <- -1*pca_clr$x[, "PC1"]
    }
    if (pca_clr$rotation[1, "PC2"] < 0) {
      pca_clr$rotation[, "PC2"] <- -1*pca_clr$rotation[, "PC2"]
      pca_clr$x[, "PC2"] <- -1*pca_clr$x[, "PC2"]
    }
    if (pca_clr$rotation[1, "PC3"] < 0) {
      pca_clr$rotation[, "PC3"] <- -1*pca_clr$rotation[, "PC3"]
      pca_clr$x[, "PC3"] <- -1*pca_clr$x[, "PC3"]
    }
    
    saveRDS(file = "output/raw_output/species_pca.Rdata")
  }
  
  
  
  
  sample_data(my_pseq) <- cbind(sample_data(my_pseq), pca_clr$x %>%
                                  as_tibble %>%
                                  select(paste0("PC", 1:n_axes)))
  
  
  return(my_pseq)
  
}

add_clr_abundances_to_meta_data <- function(pseq) {
  
  my_pseq <- pseq
  
  otus_clr <- abundances(my_pseq) %>%
    microbiome::transform("compositional") %>% 
    microbiome::transform("clr") %>% 
    t
  
  colnames(otus_clr) <- colnames(otus_clr) %>%
    gsub("\\(", "_", .) %>% 
    gsub("\\)", "", .) %>% 
    gsub(" ", "", .) 
  
  sample_data(my_pseq) <- cbind(meta(my_pseq), otus_clr)
  
  return(my_pseq)
  
  
}




add_diversities_to_meta_data <- function(pseq, species_pseq) {
  
  my_pseq <- pseq
  
  ## Diversities
  # species_pseq <- readRDS(species_level_phyloseq_path)
  species_diversities <- estimate_richness(species_pseq, measures = c("Observed", "Shannon"))
  sample_data(my_pseq) <- cbind(meta(my_pseq), species_diversities[, c("Observed", "Shannon")])
  
  return(my_pseq)
}

get_pca <- function(species_pseq, check_existence = TRUE) {
  
  set.seed(12345)
  
  
  # Check if PCA output file exists. If not, then compute
  if(file.exists("output/raw_output/species_pca.Rdata") & check_existence) {
    pca_clr <- readRDS(file = "output/raw_output/species_pca.Rdata")
  } else {
    
    print("Doing PCA. This may take a while!")
    
    pca_clr <- prcomp(microbiome::transform(
      microbiome::transform(species_pseq,
                            "compositional"),
      "clr"
    ) %>%
      abundances %>%
      t)
    
    # Flip PC1 to get a positive association with mortality
    pca_clr$rotation[, "PC1"] <- -1*pca_clr$rotation[, "PC1"]
    pca_clr$x[, "PC1"] <- -1*pca_clr$x[, "PC1"]
    
  }
  
  return(pca_clr)
  
}


# get subnet abundances, with CLR transformation clr = TRUE
subnet_abundances <- function(pseq, subnet_components, clr = "TRUE") {
  
  my_pseq <- pseq
  
  
  ## Add subnet components
  if(!exists("subnet_components")) subnet_components <- readRDS(network_identities_path)
  
  subnet_abundances <- abundances(my_pseq) %>% t %>% as.data.frame()
  
  # Aggregate subnet reads and remove the constituents
  for(i in unique(subnet_components$componentID)) {
    
    # Get subnet taxa
    subnet_taxa <- subnet_components %>%
      filter(componentID == i) %>% 
      pull(OTUunderscored)
    
    # Aggregate reads
    subnet_reads <- subnet_abundances[, subnet_taxa] %>% 
      rowSums()
    
    # Remove taxa from otu table
    subnet_abundances <- subnet_abundances[, !(colnames(subnet_abundances) %in% subnet_taxa)]
    
    # Add aggregated reads
    subnet_abundances[, paste0("subnet_", i)] <- subnet_reads
    
  }
  
  
  if(isTRUE(clr)) {
    
    # CLR transform subnet abundances
    abundances <- subnet_abundances %>%  
      t %>% 
      microbiome::transform("compositional") %>% 
      microbiome::transform("clr") %>% 
      t
    
  } else {
    abundances <- subnet_abundances
  }
  
  
  
  return(abundances)
}


# Z-normalization
z_normalize <- function(x, m = TRUE) {
  
  if(m == TRUE) {
    
    return((x-mean(x, na.rm = T))/sd(x, na.rm = T))
    
  } else {
    return(x/sd(x, na.rm = T))
    
  }
  
  
}

# Get percent of explained variation from prcomp x for axis n
PC_variation_explained <- function(x, n, round) {
  eigs <- x$sdev^2
  var <- eigs[n] / sum(eigs)
  return(100*round(var, round))
}

# wrapper: matrix to tibble with rownames to column
m_neat <- function(x, colnames = NULL) {
  x <- x %>%
    as.data.frame() %>% 
    rownames_to_column()
  
  if(!is.null(colnames)) {
    x <- x %>% set_colnames(colnames)
  }
  
  return(x)
}

# change variable names according to a df with variable explanations
change_colnames <- function(df, legend_df, col=3) {
  cols <- colnames(df)
  new_cols <- legend_df[cols, col]
  colnames(df) <- new_cols
  return(df)
}
change_varname <- function(var, legend_df, col=3) {
  row <- which(rownames(legend_df) == var)
  sober_varname <- legend_df[row, col]
  return(sober_varname)
}

# remove empty rows and columns
trim_empty_rows <- function(df, value=0) {
  
  df_reduced <- df[, colSums(df != value) > 0]
  df_reduced <- df_reduced[rowSums(df_reduced != value) > 0, ]
  
}

# clean taxa names: Bacteroides_Bacteria --> Bacteroides (Bacteria)
clean_genus_names <- function(genera, kingdom = TRUE) {
  
  if(isTRUE(kingdom)){  
    
    title <- gsub("_", " ", genera)
    title <- gsub("Bacteria", "(Bacteria)", title)
    title <- gsub("Viruses", "(Viruses)", title)
    title <- gsub("Archaea", "(Archaea)", title)
    
    for(i in 1:length(title)) {
      
      if(grepl("Plasmid", title[i])) {
        title[i] <- gsub(")Plasmid", "Plasmid)", title[i])
      }
    }
    title <- gsub("Plasmid", " Plasmid", title)
  } else {
    
    title <- gsub("_", "", genera)
    title <- gsub("Bacteria", "", title)
    title <- gsub("Viruses", " (Viruses)", title)
    title <- gsub("Archaea", " (Archaea)", title)
    title <- gsub("Plasmid", " (Plasmid)", title)
    title <- gsub("\\) \\(", " ", title)
    
    
  }
  
  
  
  return(title)
}
# dirty taxa names
dirty_genus_names <- function(x) {
  
  if(class(x) == "phyloseq") {
    genera <- taxa_names(x)
    
  } else {
    x <- genera
  }
  
  genera_edited <- genera %>% 
    gsub(" ", "_", .) %>% 
    gsub("\\(", "", .) %>% 
    sub("\\)", "", .) %>% 
    gsub("-", "_", .) %>% 
    gsub(",", "_", .)
  
  if(class(x) == "phyloseq") {
    
    taxa_names(x) <- genera_edited
    return(x)
    
  } else {
    return(genera_edited)
  }

  
}

