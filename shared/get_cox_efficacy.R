# Function to get Cox PH efficacy as in Chandramohan et al. 2021
# by year and overall 

get_cox_efficacy <- function(df, 
                             ref,
                             model = FALSE){
  
  tidy_results <- vector("list", 4)

  # Make different columns for various ref = 
  df <- df %>%
    mutate(arm_smcref = factor(arm, levels = c('smc','rtss','both','none')),
           arm_rtssref = factor(arm, levels = c('rtss','smc','both','none')),
           arm_noneref = factor(arm, levels = c('none','rtss','smc','both')),
           end_time = ifelse(start_time==end_time, end_time + 0.00001, end_time))
           
  for (i in 1:4){
    
    if(i < 4){
      d <- df %>% 
        dplyr::filter(year == i)
    } else d <- df
    
    # Create formula dynamically using the ref argument
    if(model){
      formula_str <- paste("Surv(start_time, end_time, event) ~ factor(", ref, ")")
      cox_formula <- as.formula(formula_str)
      
      clustervar = d$rid
    } else {
      formula_str <- paste("Surv(start_time, end_time, event) ~ factor(", ref, ") + factor(country)")
      cox_formula <- as.formula(formula_str)
      
      clustervar = d$rid
    }
    
    coxmodel <- coxph(cox_formula,# + factor(country), 
                                   data = d,
                                   cluster = clustervar,
                                   ties = "efron")
    
    results <- tidy(coxmodel,
                    exponentiate = TRUE,
                    conf.int = TRUE) %>%
      # Calculate vaccine efficacy and its confidence intervals
      mutate(
        term = case_when(
          term == 'factor(arm_rtssref)smc' ~ "SMC vs. RTSS",
          term == 'factor(arm_rtssref)both' ~ "Both vs. RTSS",
          term == 'factor(arm_rtssref)vaxsmc' ~ "Both vs. RTSS",
          term == 'factor(arm_rtssref)none' ~ "None vs. RTSS",
          term == 'factor(country)Mali' ~ 'Mali vs. BF',
          term == 'factor(arm_smcref)rtss' ~ "RTSS vs. SMC",
          term == 'factor(arm_smcref)both' ~ "Both vs. SMC",
          term == 'factor(arm_smcref)vax' ~ "RTSS vs. SMC",
          term == 'factor(arm_smcref)vaxsmc' ~ "Both vs. SMC",
          term == 'factor(arm_smcref)none' ~ "None vs. SMC",
          term == 'factor(arm_noneref)vax' ~ "RTSS vs None",
          term == 'factor(arm_noneref)rtss' ~ "RTSS vs None",
          term == 'factor(arm_noneref)smc' ~ "SMC vs None",
          term == 'factor(arm_noneref)vaxsmc' ~ "Both vs None",
          term == 'factor(arm_noneref)both' ~ "Both vs None",
          TRUE ~ term
        ),
        VE = 1 - estimate,              # VE = 1 - HR
        VE_lower = 1 - conf.high,       # Lower CI for VE = 1 - Upper CI for HR
        VE_upper = 1 - conf.low ,        # Upper CI for VE = 1 - Lower CI for HR
        year = ifelse(i == 4, 'overall',as.character(i))
      ) %>% mutate(
        n_events = coxmodel$nevent,
        n_obs = coxmodel$n
      )
    
    tidy_results[[i]] <- results
  }
  
  tidy_results <- bind_rows(tidy_results)
  
  return(tidy_results)
}


# to get the monthly efficacy, i think i'll need to reformat the data to be grouped by month, to have person time 
# in each month over the whole study period -- or since we are assuming no loss of follow up here, could just use 
# n as person time and have the number of cases per month / personmonths (just n)?

