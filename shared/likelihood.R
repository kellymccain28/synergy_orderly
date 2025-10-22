# Likelihood function 

calculate_poisson_likelihood <- function(trial_df = monthly_inci_trial, 
                                         model_df = monthly_inci_model,
                                         return_by_arm = TRUE){
  
  # Cases ~ Poisson(infection rate × person-months)
  # k is the observed number of cases in that month(t) and intervention group (i) (from trial)
  k <- trial_df %>%
    select(arm, yearmonth, n_cases) %>%
    pivot_wider(names_from = yearmonth, 
                values_from = n_cases,
                values_fill = 0) %>%
    column_to_rownames("arm") %>%
    as.matrix()
  
  trial_persontime <- trial_df %>%
    select(arm, yearmonth, person_months) %>%
    pivot_wider(names_from = yearmonth, 
                values_from = person_months,
                values_fill = 0) %>%
    column_to_rownames('arm') %>%
    as.matrix()
  
  
  # get modelled rates (cases per 1000 person mnths in trial)
  # lambda is the expected number of cases in that month(t) and intervention group(i) (from model) -- rate from model * person time in trial 
  model_rates <- model_df %>%
    filter(arm != 'none') %>%
    mutate(arm = case_when(
      arm == 'vaxsmc' ~ 'both',
      arm == 'vax' ~ 'rtss',
      TRUE ~ arm),
      arm = factor(arm, levels = c('both','rtss','smc'))) %>%
    select(arm, yearmonth, incidence_per_1000pm) %>%
    pivot_wider(names_from = yearmonth, 
                values_from = incidence_per_1000pm,
                values_fill = 0) %>%
    arrange(arm) %>%
    column_to_rownames("arm") %>%
    as.matrix()
  
  lambda <- model_rates * (trial_persontime / 1000)
  
  # Add tiny amount to lambda to avoid exact zeros
  epsilon <- 1e-8
  lambda <- pmax(lambda, epsilon)
  
  logL <- dpois(k, lambda, log = TRUE)
  # likelihood[i,t] <- ( lambda[i,t]^k[i,t] * exp(-lambda[i,t]) ) / factorial(k[i,t])
  # log_likelihood2 <- k * log(lambda) - lambda - lgamma(k + 1) # also works
  
  # logL <- log(likelihood)
  
  # ll <- sum(logL) # overall LL 
  
  if(return_by_arm) {
    # Return LL by arm (sum over time for each arm)
    ll_by_arm <- rowSums(logL)
    return(ll_by_arm)
  } else {
    # Return total LL (sum over arms and time)
    ll <- sum(logL)
    return(ll)
  }
  # return(ll)
}
