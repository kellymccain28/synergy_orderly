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
                               allow_superinfections = TRUE,
                               return_parasitemia = FALSE,
                               save_outputs = FALSE)
    message('finished simulation')
    
    eff <- calc_smc_efficacy(o$infection_records,
                             params_row,
                             by_week = TRUE)
    # eff_daily <- calc_smc_efficacy(o$infection_records, 
    #                                params_row, 
    #                                by_week = FALSE)
    
    eff$sim_id <- params_row$sim_id
    # eff_daily$sim_id <- params_row$sim_id
    
    matched <- observed_efficacy %>%
      left_join(eff %>% select(weeks_since_smc, efficacy) %>%
                  rename(predicted_efficacy = efficacy), by = 'weeks_since_smc')
    
    # Calculate log-likelihood
    sigma <- 0.05  # Assumed measurement error (5%)
    
    ll <- sum(dnorm(
      x = matched$observed_efficacy,          # Observed from paper
      mean = matched$predicted_efficacy,  # Predicted from model
      sd = sigma,
      log = TRUE
    ))
    
    return(-ll)
    
  },
  error = function(e) {
    warning(paste("Error in likelihood calculation:", e$message))
    return(-Inf)
  })
}
