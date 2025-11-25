calculate_efficacy_likelihood <- function(params_row,
                                          metadata_df, 
                                          base_inputs,
                                          observed_efficacy  #by week 
                                          ){
  tryCatch({
    o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                               metadata_df,
                               base_inputs,
                               output_dir = 'R:/Kelly/src/fit_smc/simulation_outputs',
                               # allow_superinfections = TRUE,
                               return_parasitemia = FALSE,
                               save_outputs = FALSE)
    message('finished simulation')
    
    eff <- calc_smc_efficacy_cumul(o$infection_records,
                                   params_row,
                                   by_week = TRUE)
    
    eff$sim_id <- params_row$sim_id
    
    matched <- observed_efficacy %>%
      mutate(weeks_since_smc = ceiling(day_since_smc / 7)) %>%
      group_by(weeks_since_smc) %>%
      summarize(observed_efficacy = mean(efficacy)) %>%
      left_join(eff %>% select(weeks_since_smc, efficacy) %>%
                  rename(predicted_efficacy = efficacy), 
                by = 'weeks_since_smc') %>% ungroup()
    
    # CHECK: Do we have matches?
    message("Matched rows: ", nrow(matched))
    
    # Remove NAs before calculating likelihood
    matched_complete <- matched %>%
      filter(!is.na(observed_efficacy), !is.na(predicted_efficacy))
    
    if(nrow(matched_complete) == 0) {
      warning("No matching timepoints between observed and predicted!")
      return(1e10)  # Return large penalty, not -Inf
    }
    
    # Calculate log-likelihood
    sigma <- 0.05  # Assumed measurement error (5%)
    
    ll <- sum(dnorm(
      x = matched_complete$observed_efficacy,          # Observed from paper
      mean = matched_complete$predicted_efficacy,  # Predicted from model
      sd = sigma,
      log = TRUE
    ))
    
    message("Log-likelihood: ", ll, " | Negative LL: ", -ll)
    
    # Check for invalid values
    if(is.na(ll) || is.infinite(ll)) {
      warning("Invalid log-likelihood calculated: ", ll)
      return(1e10)
    }
    
    return(-ll) # we will minimize the negative log likelihood
    
  },
  error = function(e) {
    warning(paste("Error in likelihood calculation:", e$message))
    return(1e10)
  })
}


calculate_efficacy_likelihood_rtss <- function(params_row,
                                               metadata_df, 
                                               base_inputs,
                                               observed_efficacy_rtss  #by week 
){
  tryCatch({
    o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                               metadata_df,
                               base_inputs,
                               output_dir = 'R:/Kelly/src/fit_smc/simulation_outputs',
                               # allow_superinfections = TRUE,
                               return_parasitemia = FALSE,
                               save_outputs = FALSE)
    message('finished simulation')
    
    o$infection_records$sim_id <- params_row$sim_id
    
    # Follow up begins 3 weeks after 3rd dose in trial to which White 2015 fit 
    # this effectively removes any cases among the RTSS arm that would have 
    # been bitten/infected before follow up 
    # so to approximate this, will remove infectiosn that happened before follow-up
    infs <- o$infection_records %>%
      filter(BSinfection_day >= 0)
    
    eff <- calc_rtss_efficacy(infs)
    
    matched <- observed_efficacy_rtss %>%
      left_join(eff %>% select(weeks_since_rtss, efficacy) %>%
                  rename(predicted_efficacy = efficacy), 
                by = 'weeks_since_rtss')
    
    # CHECK: Do we have matches?
    message("Matched rows: ", nrow(matched))
    
    # Remove NAs before calculating likelihood
    matched_complete <- matched %>%
      filter(!is.na(observed_efficacy), !is.na(predicted_efficacy))
    
    if(nrow(matched_complete) == 0) {
      warning("No matching timepoints between observed and predicted!")
      return(1e10)  # Return large penalty, not -Inf
    }
    
    # Calculate log-likelihood
    sigma <- 0.05  # Assumed measurement error (5%)
    
    ll <- sum(dnorm(
      x = matched_complete$observed_efficacy,          # Observed from paper
      mean = matched_complete$predicted_efficacy,  # Predicted from model
      sd = sigma,
      log = TRUE
    ))
    
    message("Log-likelihood: ", ll, " | Negative LL: ", -ll)
    
    # Check for invalid values
    if(is.na(ll) || is.infinite(ll)) {
      warning("Invalid log-likelihood calculated: ", ll)
      return(1e10)
    }
    
    return(-ll) # we will minimize the negative log likelihood
    
  },
  error = function(e) {
    warning(paste("Error in likelihood calculation:", e$message, e$warning))
    return(1e10)
  })
}

