calculate_efficacy_likelihood <- function(params_row,
                                          metadata_df, 
                                          base_inputs,
                                          observed_efficacy  #by week 
                                          ){
  tryCatch({
    startsim <- Sys.time()
    o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                               metadata_df,
                               base_inputs,
                               output_dir = 'R:/Kelly/src/fit_smc/simulation_outputs',
                               # allow_superinfections = TRUE,
                               return_parasitemia = TRUE,
                               save_outputs = FALSE)
    endsim <- Sys.time()
    message('time to do sim: ', endsim - startsim)
    message('finished simulation')
    filt <- o$infection_records %>% filter(time_ext >=0) # filtering out the bites that happened before the start of FU
    eff <- calc_smc_efficacy_cumul(filt,#o$infection_records,
                                   params_row,
                                   by_week = TRUE)
    
    eff$sim_id <- params_row$sim_id
    
    matched <- observed_efficacy %>%
      mutate(weeks_since_smc = ceiling(day_since_smc / 7)) %>%
      group_by(weeks_since_smc) %>%
      summarize(observed_efficacy = mean(efficacy)) %>%
      mutate(observed_efficacy = ifelse(observed_efficacy == 1, 0.999, observed_efficacy)) %>%
      left_join(eff %>% select(weeks_since_smc, efficacy) %>%
                  rename(predicted_efficacy = efficacy), 
                by = 'weeks_since_smc') %>% ungroup()
    
    # CHECK: Do we have matches?
    message("Matched rows: ", nrow(matched))
    
    # Remove NAs before calculating likelihood
    matched_complete <- matched %>%
      filter(!is.na(observed_efficacy), !is.na(predicted_efficacy))
    
    # if(nrow(matched_complete) == 0) {
    #   warning("No matching timepoints between observed and predicted!")
    #   return(1e10)  # Return large penalty, not -Inf
    # }
    # 
    # Calculate mean least squares
    mls <- mean((matched_complete$observed_efficacy - matched_complete$predicted_efficacy)^2)
    
    message("Mean Least Squares: ", mls)
    
    # # Calculate log-likelihood
    # sigma <- 0.05  # Assumed measurement error (5%)
    # 
    # logit_obs <- qlogis(matched_complete$observed_efficacy)
    # logit_pred <- qlogis(matched_complete$predicted_efficacy)
    # 
    # ll <- sum(dnorm(
    #   x = logit_obs[2:9],#matched_complete$observed_efficacy,          # Observed from paper
    #   mean = logit_pred[2:9],#matched_complete$predicted_efficacy,  # Predicted from model
    #   sd = sigma,
    #   log = TRUE
    # ))
    # 
    # message("Log-likelihood: ", ll, " | Negative LL: ", -ll)
    # 
    # Check for invalid values
    if(is.na(mls) || is.infinite(mls)) {
      warning("Invalid MLS calculated: ", mls)
      return(1e10)
    }
    # 
    # return(-ll) # we will minimize the negative log likelihood
    
    return(mls)
  },
  error = function(e) {
    warning(paste("Error:", e$message))
    # return(1e10)
  })
}


calculate_efficacy_likelihood_rtss <- function(params_row,
                                               metadata_df, 
                                               base_inputs,
                                               observed_efficacy_rtss  #by week 
){
  tryCatch({
    message('starting simulation now')
    
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
    
    # if(nrow(matched_complete) == 0) {
    #   warning("No matching timepoints between observed and predicted!")
    #   return(1e10)  # Return large penalty, not -Inf
    # }
    
    # # Calculate log-likelihood
    # sigma <- 0.05  # Assumed measurement error (5%)
    # 
    # ll <- sum(dnorm(
    #   x = matched_complete$observed_efficacy,          # Observed from paper
    #   mean = matched_complete$predicted_efficacy,  # Predicted from model
    #   sd = sigma,
    #   log = TRUE
    # ))
    # 
    # message("Log-likelihood: ", ll, " | Negative LL: ", -ll)
    # 
    
    # 
    # return(-ll) # we will minimize the negative log likelihood
    
    # Calculate mean least squares
    mls <- mean((matched_complete$observed_efficacy - matched_complete$predicted_efficacy)^2)
    
    message("Mean Least Squares: ", mls)
    
    # # Check for invalid values
    if(is.na(mls) || is.infinite(mls)) {
      warning("Invalid MLS calculated: ", mls)
      return(1e10)
    }
    return(mls)
    
  },
  error = function(e) {
    warning(paste("Error:", e$message, e$warning))
    # return(1e10)
  })
}

