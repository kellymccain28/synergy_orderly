get_cox_efficacy <- function(df, 
                             ref,
                             model = FALSE){
  
  tidy_results <- vector("list", 4)
  
  # Create formula dynamically using the ref argument
  if(model){
    formula_str <- paste("Surv(start_time, end_time, event) ~ factor(", ref, ")")
    cox_formula <- as.formula(formula_str)
  } else {
    formula_str <- paste("Surv(start_time, end_time, event) ~ factor(", ref, ") + factor(country)")
    cox_formula <- as.formula(formula_str)
  }
  
  for (i in 1:4){
    
    if(i < 4){
      d <- df %>% 
        filter(year == i)
    } else d <- df
    
    coxmodel <- coxph(cox_formula,# + factor(country), 
                                   data = d,
                                   cluster = rid,
                                   ties = "efron")
    
    results <- tidy(coxmodel,
                    exponentiate = TRUE,
                    conf.int = TRUE) %>%
      # Calculate vaccine efficacy and its confidence intervals
      mutate(
        term = case_when(
          term == 'factor(arm_rtssref)smc' ~ "SMC vs. RTSS",
          term == 'factor(arm_rtssref)both' ~ "Both vs. RTSS",
          term == 'factor(country)Mali' ~ 'Mali vs. BF',
          term == 'factor(arm_smcref)vax' ~ "RTSS vs. SMC",
          term == 'factor(arm_smcref)vaxsmc' ~ "Both vs. SMC",
          term == 'factor(arm_smcref)none' ~ "None vs. SMC",
          TRUE ~ term
        ),
        VE = 1 - estimate,              # VE = 1 - HR
        VE_lower = 1 - conf.high,       # Lower CI for VE = 1 - Upper CI for HR
        VE_upper = 1 - conf.low ,        # Upper CI for VE = 1 - Lower CI for HR
        year = ifelse(i == 4, 'overall',as.character(i))
      ) %>% mutate(
        n_events = ag_model_rtssref_year$nevent,
        n_obs = ag_model_rtssref_year$n
      )
    
    tidy_results[[i]] <- results
  }
  
  tidy_results <- bind_rows(tidy_results)
  
  return(tidy_results)
}
