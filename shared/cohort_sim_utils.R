# functions to do parameter sweeps 

# Main cohort simulation function
run_cohort_simulation <- function(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                                  metadata_df, 
                                  base_inputs, 
                                  return_parasitemia = TRUE,
                                  output_dir = "R:/Kelly/src/fit_smc/simulation_outputs",
                                  save_outputs = FALSE) {
  
  message('running cohort simulation for: ', params_row$sim_id)
  # Extract parameters from the row
  max_SMC_kill_rate <- params_row$max_SMC_kill_rate
  smc_lambda <- params_row$lambda #17
  smc_kappa <- params_row$kappa #0.28 - these are default values (can't rmember where I got them though..)
  p_bite <- unlist(params_row$p_bite)
  alpha_ab_value <- if(any(names(params_row) == 'alpha_ab')) params_row$alpha_ab else 1.32
  beta_ab_value <- if(any(names(params_row) == 'beta_ab')) params_row$beta_ab else 6.62
  vmin_value <- if(any(names(params_row) == 'vmin')) params_row$vmin else 0
  
  message("Running simulation: ", params_row$sim_id)
  
  # Extract base inputs (these stay constant across parameter sets)
  trial_ts <- base_inputs$trial_timesteps#[[i]]
  burnin <- base_inputs$burnin
  threshold <- base_inputs$threshold
  VB <- base_inputs$VB
  tstep <- base_inputs$tstep
  t_liverstage <- base_inputs$t_liverstage
  country_to_run <- base_inputs$country
  country_short <- base_inputs$country_short
  
  # Used when allow_superinfections = FALSE to move treated people to be no longer immune
  treatment_probability <- base_inputs$treatment_probability
  successful_treatment_probability <- base_inputs$successful_treatment_probability
  
  # now, just do 100 timesteps (200 days) for each infection, which is plenty of time 
  ts <- trial_ts
  tt <- seq(1, ts, by = tstep)#seq(0, ((trial_ts + burnin) - t)/divide, by = tstep)
  
  # Initialize storage for this simulation
  infection_records <- data.frame(
    rid = integer(),
    BSinfection_day = integer(),
    threshold_day = integer(),
    t_toreach_threshold = integer(),
    vaccination_day = integer(),
    PEV = integer(),
    SMC = integer(),
    arm = character(),
    p_bite = double(),
    recovery_day = double(),
    country = character(),
    receives_treatment = logical(),
    treatment_efficacy = logical(),
    treatment_successful = logical(), 
    treatment_day = integer()
  )
  
  N <- nrow(metadata_df)
  
  # Make weights to prevent homogeneous transmission 
  weights <- rexp(N) # Generate random weights to sample children at different probabilities (could also use rpois(N) or rexp(N))
  weights <- weights / sum(weights) # normalize to sum to 1
  
  message('calculating time since smc and parasite kill rate vector')
  # Calculate SMC kill rate vectors for each child  - not sure where these are from 
  # this is time since smc where the day is relative to the first day of the simulation (not including burnin)
  if("smc_dose_days" %in% names(params_row)){ 
    metadata_df$smc_dose_days <- params_row$smc_dose_days
  } 
  
  metadata_df$smckillvec <- lapply(metadata_df$smc_dose_days, get_smc_vectors,
                                   ts = trial_ts,
                                   burnin = base_inputs$burnin,
                                   max_SMC_kill_rate = max_SMC_kill_rate,
                                   lambda = smc_lambda,
                                   kappa = smc_kappa)
  
  # hill? plot(1 - (1/(1+(3.7/seq(0,6,0.1))^4))) plot(1 - (1/(1+(25/seq(0,60,1))^3.2)))
  message('Calculated smc kill rate vectors for each child')
  # this kill vector has 25 2-day timesteps for burnin=50, then 548 2-day timesteps for the 1095 day simulation  
  
  # Initialize parasitemia storage
  parasitemia_storage <- vector("list", length = trial_ts + burnin)
  
  susceptibles <- rep(TRUE, N)
  
  message('starting main sim loop')
  # Main simulation loop 
  for (t in 1:(trial_ts + burnin)) {
    # Determine number of new infectious bites on time t
    n_infectious_bites <- rpois(1, lambda = N * p_bite[t])
    # message("probability of a bite:", p_bite[t])
    if (n_infectious_bites == 0) next
    
    # Sample random children to be bitten 
    bites <- sample(metadata_df$rid, size = n_infectious_bites, prob = weights, replace = TRUE)
    
    bit_kids <- unique(bites)
    # Remove any bit kids with no vaccine follow-up data 
    bit_kids <- bit_kids[bit_kids %in% metadata_df$rid] 
    
    # Update susceptibility vector so that recovered kids would be susceptible again
    if (exists("infection_records") && nrow(infection_records) > 0 ){#& !allow_superinfections) {
      
      # Filter to only records with actual detections
      detected_records <- infection_records[!is.na(infection_records$detection_day), ]
      
      if (nrow(detected_records) > 0) {
        # Get the most recent detection for each person
        most_recent <- detected_records %>%
          group_by(rid) %>%
          slice_max(detection_day, n = 1) %>%  # Get most recent detection_day for each rid
          ungroup()
        
        # Find who should recover today (day t)
        recovery_today <- (most_recent$recovery_day + burnin) == t #????? & most_recent$treatment_successful???? #add burnin because the recovery day calculation is the external time but for this siulation we want ot know if on time t they could be infected
        # if the recovery day equals today then these are the ones that are newly susceptible 
        if (any(recovery_today)) {
          recovering_kids <- most_recent$rid[recovery_today]
          susceptibles[recovering_kids] <- TRUE
        }
        # message(cat(c('recovered today: ', most_recent$rid[recovery_today]), sep = ','))
        
        # Successfully treated children are not susceptible from treatment day until recovery
        treated_and_protected <- most_recent$treatment_successful & 
          (t >= (most_recent$treatment_day + burnin)) & 
          (t < (most_recent$recovery_day + burnin))
        if (any(treated_and_protected)) {
          protected_kids <- most_recent$rid[treated_and_protected]
          susceptibles[protected_kids] <- FALSE
        }
        # message(cat("treated and protected: ", most_recent$rid[treated_and_protected], sep = ','))
      }
    }
    
    # Avoid superinfection -- if a kid is bitten and they aren't susceptible they are not added to the list of bit_kids
    # bit_kids_t <- c()
    # j <- 1
    # for (kid in bit_kids){
    #   if (susceptibles[kid]){ 
    #     bit_kids_t[j] <- kid # Add this kid to the filtered list if the kid is susceptible
    #     
    #     j <- j + 1 # move to the next position in the filtered list 
    #   }
    # }
    bit_kids_susceptible <- bit_kids[susceptibles[bit_kids]]
    
    if(t %% 10 == 0){
      message(str_glue('{length(bit_kids_susceptible)} / {length(bit_kids)} bit kids were susceptible on time ', t, ' out of ', (trial_ts+burnin)))
    }
    bit_kids <- bit_kids_susceptible
    rm(bit_kids_susceptible)
    
    if(length(bit_kids) > 0) {
      kid_metadata <- metadata_df[metadata_df$rid %in% bit_kids, ]
      
      if(t < burnin) {
        # No interventions during burnin
        PEV_vec <- rep(0, length(bit_kids))
        # these both need a placeholder value (won't matter if PEV = 0 or SMC = 0)
        t_since_vax_vec <- rep(0, length(bit_kids))
        t_toboost1_vec <- rep(500, length(bit_kids))
        t_toboost2_vec <- rep(1000, length(bit_kids))
        
        # Because the SMC vector is already made, with a padded start for the burnin, we need to have SMC on in the model for the SMC children 
        # Make named vector of SMC and map to bit_kids order 
        smc_lookup <- setNames(kid_metadata$SMC, kid_metadata$rid)
        SMC_vec <- smc_lookup[as.character(bit_kids)]
        
        # Get the kill rate vectors for SMC
        # Find the row index for the target rid
        smc_kill_vec_lookup <- setNames(kid_metadata$smckillvec, kid_metadata$rid)
        SMC_kill_vec <- smc_kill_vec_lookup[as.character(bit_kids)]
        # subset the kill rate vector to be from the external time (t includes the burnin) to the end of the vector
        # vector is every two days 
        SMC_kill_vec <- lapply(SMC_kill_vec, function(smcvec){
          subset <- smcvec[floor((t + t_liverstage) / 2) :length(smcvec)] # t_liverstage added because we want to subset from when BS infection starts + the liver stage time, to the end; +burnin removed because it is already taken into account in the t
          return(subset)
        })
        
        SMC_timev <- rep(list(0:(length(SMC_kill_vec[[1]])-1)), 
                         length = length(bit_kids))
      } else {
        
        # Create named vectors of PEVand vaccination days 
        pev_lookup <- setNames(kid_metadata$PEV, kid_metadata$rid)
        vax_day_lookup <- setNames(kid_metadata$vaccination_day, kid_metadata$rid)
        t_toboost1_lookup <- setNames(kid_metadata$t_to_boost1, kid_metadata$rid)
        t_toboost2_lookup <- setNames(kid_metadata$t_to_boost2, kid_metadata$rid)
        
        # Map to bit_kids order
        PEV_vec <- pev_lookup[as.character(bit_kids)]
        # Find the time since vaccination -- time between 3rd dose (vax day) and infectious bite
        t_since_vax_vec <- (t - burnin) - vax_day_lookup[as.character(bit_kids)] # a (-) value of vaxdaylookup indicates that vaccination was after BS infection; t is infectious bite 
        # vaccination before the infection begins 
        # for negative values of t_since_vax_vec, which means that the vaccine was delivered after the current time t, change to 0 so that ab_user is not NA
        t_since_vax_vec <- ifelse(t_since_vax_vec<0, 0, t_since_vax_vec)
        
        # booster doses -- these are the timings of boosters 1 and 2 relative to the 3rd dose - used in ab calculation
        t_toboost1_vec <- t_toboost1_lookup[as.character(bit_kids)]
        t_toboost2_vec <- t_toboost2_lookup[as.character(bit_kids)]
        
        # Make named vector of SMC and map to bit_kids order 
        smc_lookup <- setNames(kid_metadata$SMC, kid_metadata$rid)
        SMC_vec <- smc_lookup[as.character(bit_kids)]
        
        # Get the kill rate vectors for SMC 
        # Find the row index for the target rid
        smc_kill_vec_lookup <- setNames(kid_metadata$smckillvec, kid_metadata$rid)
        SMC_kill_vec <- smc_kill_vec_lookup[as.character(bit_kids)]
        # subset the kill rate vector to be from the external time (t includes the burnin) to the end of the vector
        # vector is every two days 
        SMC_kill_vec <- lapply(SMC_kill_vec, function(smcvec){
          subset <- smcvec[floor((t + t_liverstage) / 2) :length(smcvec)] # t_liverstage added because we want to subset from when BS infection starts which is the bite + the liver stage time, to the end; +burnin removed because it is already taken into account in the t
          return(subset)
        })
        
        SMC_timev <- rep(list(0:(length(SMC_kill_vec[[1]])-1)), 
                         length = length(bit_kids))
        
      }
      
      # Find time length until end of follow-up to run infection sim  
      tt_until_end_cohort <- seq(1, length(SMC_kill_vec[[1]]))#seq(1, ((trial_ts + burnin) - t)/2, by = tstep)#seq(1, 100, by = tstep)#
      # message("timesteps to run individual sims:", tt_until_end_cohort)
      # message('smc vector: ', SMC_kill_vec[[1]])
      
      # Run within-host simulation with current parameters
      params_df <- data.frame(
        PEV_on = PEV_vec,
        SMC_on = SMC_vec,
        t_inf_vax = t_since_vax_vec,
        infection_start_day = t - burnin + t_liverstage, # (day of bite is t-burnin and infectionstartdayis BS infection start day) because in run_process_model it is only used for formatting the time var -- time*2 + infstartday so that value must be without the burnin period 
        SMC_time = I(SMC_timev),
        SMC_kill_vec = I(SMC_kill_vec), 
        rid = bit_kids,
        alpha_ab = alpha_ab_value, # default values from White 2013 
        beta_ab = beta_ab_value, # default values from White 2013
        vmin = vmin_value,
        t_toboost1 = t_toboost1_vec,
        t_toboost2 = t_toboost2_vec
      )
      
      outputs <- pmap(params_df, 
                      function(PEV_on, 
                               SMC_on, 
                               t_inf_vax,  # time since vaccination -- influences ab titre 
                               infection_start_day, # current cohort day -- influences the formatting -- time variable 
                               SMC_time,
                               SMC_kill_vec,
                               rid,
                               alpha_ab, 
                               beta_ab,
                               vmin,
                               t_toboost1,
                               t_toboost2) {
                        result <- run_process_model(
                          PEV_on = PEV_on,
                          SMC_on = SMC_on,
                          t_inf_vax = t_inf_vax,
                          VB = VB,
                          tt = tt_until_end_cohort, # run the model until the end of the follow-up
                          infection_start_day = infection_start_day,# external time that infection begins 
                          SMC_time = unlist(SMC_time),
                          SMC_kill_vec = unlist(SMC_kill_vec),
                          alpha_ab = alpha_ab, # default values from White 2013 
                          beta_ab = beta_ab, # default values from White 2013
                          vmin = vmin, 
                          tboost1 = t_toboost1,
                          tboost2 = t_toboost2
                        )
                        
                        result$trajectory$rid <- rid
                        
                        return(result)
                      })
      
      # Vectorized data frame creation of infection records
      # new_records <- data.frame(
      #   rid = bit_kids,
      #   time_ext = rep(t - burnin, length(bit_kids)),                                                      # external cohort time (t is the cohort time, then we wnat to scale to be + if after burnin)
      #   t = t,  # simulation day
      #   infectious_bite_day = rep((t - burnin) - t_liverstage, length(bit_kids)),                          # bitten on day t, then assuming the liver stage takes t_liverstage days
      #   BSinfection_day = rep((t - burnin), length(bit_kids)),                                             # after liver stage, the BS begins #+ t_liverstage
      #   threshold_day = sapply(outputs, function(x) x$threshold_day),               # days since BS starts that threshold is reached
      #   detection_day = (t - burnin) + sapply(outputs, function(x) x$threshold_day),# day in cohort simulation that threshold is reached, threshold is the day since BS infection that reaches threshold #+ t_liverstage 
      #   t_toreach_threshold = sapply(outputs, function(x) x$threshold_day) + t_liverstage,   # time to reach threshold value / detection since the bite  
      #   vaccination_day =  kid_metadata$vaccination_day, #if(t < burnin) rep(NA, length(bit_kids)) else  # day of vaccination relative to the start of follow-up (day 0 external time)
      #   prob_bite = rep(p_bite[t], length(bit_kids)),
      #   recovery_day = ((t - burnin) + sapply(outputs, function(x) x$threshold_day)) + 12 - t_liverstage, # day that the child would be 'recovered' if we assume that a child is treated and has a period of prophylaxis for 12 days after detection day 
      #   # (90% at 12 days in paper but here, assuming 100% for 12 days) after the day of treatment and that all infectiosn are treated with AL 10.1038/ncomms6606
      #   country = country_to_run
      # ) 
      # Vectorized data frame creation of infection records with infectious bite day as central day 
      new_records <- data.frame(
        rid = bit_kids,
        time_ext = rep(t - burnin, length(bit_kids)),                                                      # external cohort time (t is the cohort time, then we wnat to scale to be + if after burnin)
        t = t,  # simulation day
        infectious_bite_day = rep((t - burnin), length(bit_kids)),                          # bitten on day t, then assuming the liver stage takes t_liverstage days
        BSinfection_day = rep((t - burnin) + t_liverstage, length(bit_kids)),                                             # after liver stage, the BS begins #+ t_liverstage
        threshold_day = sapply(outputs, function(x) x$threshold_day),               # days since BS starts that threshold is reached
        detection_day = (t - burnin) + sapply(outputs, function(x) x$threshold_day) + t_liverstage,# day in cohort simulation that threshold is reached, threshold is the day since BS infection that reaches threshold #+ t_liverstage 
        t_toreach_threshold = sapply(outputs, function(x) x$threshold_day) + t_liverstage,   # time to reach threshold value / detection since the bite  
        vaccination_day =  kid_metadata$vaccination_day, #if(t < burnin) rep(NA, length(bit_kids)) else  # day of vaccination relative to the start of follow-up (day 0 external time)
        prob_bite = rep(p_bite[t], length(bit_kids)),
        recovery_day = ((t - burnin) + sapply(outputs, function(x) x$threshold_day)) + 12, # day that the child would be 'recovered' if we assume that a child is treated and has a period of prophylaxis for 12 days after detection day 
        # (90% at 12 days in paper but here, assuming 100% for 12 days) after the day of treatment and that all infectiosn are treated with AL 10.1038/ncomms6606 (this is detection day  + 12 - liverstage, but detection day is threshold + t + liverstage, so the lviverstages cancel out leaving us with 12 +threshold + t)
        country = country_to_run
      ) 
      
      detected_kids <- !is.na(new_records$detection_day)#bit_kids[!is.na(new_records$detection_day)]
      
      # Initialize treatment columns (all FALSE/NA for non-detected)
      receives_treatment <- rep(FALSE, length(bit_kids))
      treatment_efficacy <- rep(FALSE, length(bit_kids))
      treatment_successful <- rep(FALSE, length(bit_kids))
      treatment_day <- rep(NA, length(bit_kids))
      
     if (any(detected_kids)) {
        receives_treatment[detected_kids] <- runif(sum(detected_kids)) < treatment_probability
        treatment_efficacy[detected_kids] <- runif(sum(detected_kids)) < successful_treatment_probability
        treatment_successful[detected_kids] <- receives_treatment[detected_kids] & treatment_efficacy[detected_kids]
        treatment_day[detected_kids] <- new_records$detection_day[detected_kids]
     }
      
      # Add treatment columns to new_records
      new_records$receives_treatment <- receives_treatment
      new_records$treatment_efficacy <- treatment_efficacy
      new_records$treatment_successful <- treatment_successful
      new_records$treatment_day <- ifelse(new_records$treatment_successful, treatment_day, NA) # treatment day is detection day -- which is t - burnin +thresholdday
      
      # Now, remove any infections that are  within 7 days of the previous infeciton 
      # find the new detectable rids 
      # new_rids <- new_records[!is.na(new_records$detection_day),] %>% arrange(rid)
      # existing_infections <- infection_records[infection_records$rid %in% new_rids$rid,]
      # if(nrow(existing_infections) > 0){
      #   most_recent_old_infections <- existing_infections %>%
      #     group_by(rid) %>%
      #     slice_max(detection_day, n = 1) %>%  # Get most recent detection_day for each rid
      #     ungroup()
      #   
      #   new_rids$toosoon <- new_rids$detection_day - most_recent_old_infections$detection_day >= 7
      #   
      #   #filter out those that are too soon (<7 days since last det day)
      #   new <- new_records %>% filter(!(rid %in% new_rids[new_rids$toosoon == FALSE]))
      #   
      #   # check that the infections in the most recent old ones are at least 7 days before the new ones in this timestep 
      #   
      # }
      
      infection_records <- rbind(infection_records, new_records)
      
      # Only update the susceptibility vector if allow_superinfections == FALSE
      # Mark as non-susceptible only those who were detected AND successfully treated
      kids_treated_successfully <- bit_kids[treatment_successful]
      
      if (length(kids_treated_successfully) > 0) {
        susceptibles[kids_treated_successfully] <- FALSE
      } 
      
      
      if(return_parasitemia){
        # Vectorized parasitemia storage creation
        parasitemia_data <- map2(outputs, bit_kids, function(output, kid) {
          
          output$trajectory %>%
            mutate(
              day1_BSinfection = t - burnin + t_liverstage,
              threshold_day = output$threshold_day,
              detection_day = t -burnin + threshold_day + t_liverstage,# do not need to multiply threshold day by 2 here since already done in run_process_model #+ t_liverstage 
              # time_ext = time_orig*2 + (t - burnin) - 1,# + infection_start_day/2,# external time should be dependent on when infection was, relative to external time(inf_start_day aka t); time (model time) has already been multiplied by 2 and related to infection time  #time*2 + (t - 1) - burnin,#+ t_liverstage
              arm = metadata_df[metadata_df$rid == kid, ]$arm,
              t = t
            )
        })
        
        
        # Store all parasitemia data for this time step
        parasitemia_storage[[t]] <- parasitemia_data
      }
    }
  }
  
  message("finished sim")
  
  if(return_parasitemia){
    # Process final outputs
    parasitemia_storage <- unlist(parasitemia_storage, recursive = FALSE)
    
    parasitemia_df <- bind_rows(parasitemia_storage) %>%
      # don't keep any infections from before follow-up begins
      filter(detection_day > 0 | is.na(detection_day)) %>%
      group_by(rid, day1_BSinfection) %>%
      mutate(
        #   child_dayinf = paste0(rid, ", day ", day1_BSinfection),
        det = ifelse(!is.na(detection_day), 1, 0),
        country = country_to_run
      ) 
    
    # keep a small sample of the parasitemia dataset
    parasitemia_df <- parasitemia_df[parasitemia_df$rid %in%
                                       (parasitemia_df %>%
                                          distinct(rid, arm) %>%
                                          group_by(arm) %>%
                                          slice_sample(n=30) %>%
                                          pull(rid)), ]
  }
  
  # Process infection records
  child_counts <- metadata_df %>%
    distinct(rid, arm, country) %>%
    count(arm, country, name = "children_in_group")
  
  infection_records_final <- infection_records %>%
    mutate(rid_original = paste0(country_short, sprintf("%04d", rid))) %>%
    filter(detection_day > 0 | is.na(detection_day)) %>% # no detection days before follow-up starts 
    left_join(metadata_df %>% select(-vaccination_day), 
              by = c('rid_original', 'rid', 'country')) %>%
    left_join(child_counts, by = c('arm', 'country')) %>%
    mutate(detectable = ifelse(is.na(threshold_day), 0, 1)) %>%
    group_by(arm) %>%
    arrange(detection_day) %>%
    mutate(cumul_inf = cumsum(detectable)) %>%
    filter(as.Date(detection_day, origin = as.Date('2017-04-01')) > v1_date | is.na(detection_day)) # remove any infections that started before v1_date (beginning of follow-up for each individual in trial)
  
  message('processed sim')
  
  # Save outputs if requested
  if (save_outputs && !is.null(output_dir)) {
    # Create diagnostic plots
    plots <- create_diagnostic_plots(infection_records_final,
                                     max_SMC_kill_rate,
                                     smc_lambda, smc_kappa)
    message('made plots')
    
    if(return_parasitemia){
      # Prepare results
      results <- list(
        sim_id = params_row$sim_id,
        parameters = params_row,
        infection_records = infection_records_final,
        parasitemia_data = parasitemia_df,
        diagnostic_plots = plots
      )
    } else {
      # Prepare results
      results <- list(
        sim_id = params_row$sim_id,
        parameters = params_row,
        infection_records = infection_records_final,
        diagnostic_plots = plots)
    }
    
    save_simulation_outputs(results, output_dir)
    message('saved outputs')
    rm(results)
  }
  
  if(!save_outputs){
    if(return_parasitemia){
      # Prepare results
      results <- list(
        sim_id = params_row$sim_id,
        parameters = params_row,
        infection_records = infection_records_final,
        parasitemia_data = parasitemia_df
      )
    } else {
      results <- list(
        sim_id = params_row$sim_id,
        parameters = params_row,
        infection_records = infection_records_final
      )
    }
    
    return(results)
  }
}


# Function to create diagnostic plots
create_diagnostic_plots <- function(#parasitemia_df, 
  infection_records_final,
  max_SMC_kill_rate, smc_lambda, smc_kappa) {
  
  # Plot 1: Parasitemia over time
  # p1 <- ggplot(parasitemia_df) +
  #   geom_line(aes(x=time_ext/365.25, y = parasites, 
  #                 group = as.factor(paste0(rid, ", day ", day1_BSinfection)), 
  #                 color = det), alpha = 0.5) +
  #   scale_y_log10() +
  #   geom_hline(aes(yintercept = 10), linetype = 2, color = 'darkred', linewidth = 1) +
  #   geom_hline(aes(yintercept = 1e-5), linetype = 2, color = 'darkgreen', linewidth = 1) +
  #   geom_vline(aes(xintercept = 0), color = 'darkturquoise', linetype = 2, linewidth = 1) +
  #   facet_wrap(~arm) +
  #   theme_bw() +
  #   theme(legend.position = 'none') +
  #   labs(title = "Parasitemia trajectories",
  #        caption = str_glue("SMC max kill rate: {max_SMC_kill_rate}, SMC lambda: {smc_lambda}, SMC kappa: {smc_kappa}"))
  
  # Plot 2: Proportion detectable
  p2 <- ggplot(infection_records_final) +
    geom_bar(aes(x = arm, group = as.factor(detectable), 
                 fill = as.factor(detectable)), position = 'fill') +
    scale_fill_manual(values = c('darkmagenta','goldenrod')) +
    labs(fill = 'Detectable', title = "Proportion detectable infections",
         caption = str_glue("SMC max kill rate: {max_SMC_kill_rate}, SMC lambda: {smc_lambda}, SMC kappa: {smc_kappa}"))
  
  # Plot 3: Incidence
  p3 <- ggplot(infection_records_final) +
    geom_histogram(aes(x = detection_day, color=arm), stat = 'count') +
    facet_wrap(~arm) +
    labs(y = 'N infections', x = 'Day since start of follow up period',
         title = "Incidence over time") +
    theme(legend.position = 'none',
          caption = str_glue("SMC max kill rate: {max_SMC_kill_rate}, SMC lambda: {smc_lambda}, SMC kappa: {smc_kappa}"))
  
  # Plot 4: Cumulative infections
  p4 <- ggplot(infection_records_final) +
    geom_line(aes(x = detection_day, y = cumul_inf, color=arm)) +
    labs(y = 'N infections', 
         x = 'Day since start of follow up period',
         title = "Cumulative infections",
         caption = str_glue("SMC max kill rate: {max_SMC_kill_rate}, SMC lambda: {smc_lambda}, SMC kappa: {smc_kappa}"))
  
  p5 <- ggplot(infection_records_final) + 
    geom_line(aes(x = detection_day, y = prob_bite, color= arm)) +
    labs(y = 'Probability of an infectious bite', 
         x = 'Day since start of follow up period',
         title = "Probability of an infectious bite")
  
  return(list(
    # parasitemia_trajectories = p1,
    proportion_detectable = p2,
    incidence = p3,
    cumulative_infections = p4,
    prob_infectious_bite = p5
  ))
}

# Function to save simulation outputs
save_simulation_outputs <- function(results, output_dir) {
  
  # dir.create('outputs/', showWarnings = FALSE)
  
  # sim_dir <- file.path(output_dir, results$sim_id)
  
  # if (!dir.exists(sim_dir)) {
  #   dir.create(sim_dir, recursive = TRUE)
  # }
  if (!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }
  
  plots_dir <- file.path("/plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir)
  }
  
  # Save input files 
  # write.csv(params_df, file.path(sim_dir, "parameter_grid.csv"), row.names = FALSE)
  saveRDS(base_inputs, file.path(output_dir, "/base_inputs.rds"))
  
  # Save data files
  saveRDS(results$infection_records, paste0(output_dir, "/infection_records",results$sim_id,".rds"))
  if('parasitemia_data' %in% names(results)){
    saveRDS(results$parasitemia_data, paste0(output_dir, "/parasitemia_data",results$sim_id,".rds"))
  } 
  saveRDS(results$parameters %>% select(-p_bite), paste0(output_dir, "/parameters", results$sim_id, ".rds"))
  
  # Save plots
  ggsave(paste0(plots_dir, "proportion_detectable",results$sim_id,".png"), 
         results$diagnostic_plots$proportion_detectable, 
         height = 8, width = 12)
  ggsave(paste0(plots_dir, "/incidence",results$sim_id,".png"), 
         results$diagnostic_plots$incidence, 
         height = 6, width = 14)
  ggsave(paste0(plots_dir, "/cumulative_infections",results$sim_id,".png"), 
         results$diagnostic_plots$cumulative_infections, 
         height = 8, width = 12)
  ggsave(paste0(plots_dir, "/prob_infectious_bite",results$sim_id,".png"), 
         results$diagnostic_plots$prob_infectious_bite, 
         height = 8, width = 12)
  
  message("Saved outputs for ", results$sim_id)
}

# function to calculate efficacy of RTSS versus no intervention 
calc_rtss_efficacy <- function(df){
  n_rtss <- nrow(metadata_df[metadata_df$arm == 'rtss',])
  n_none <- nrow(metadata_df[metadata_df$arm == 'none',])
  vax_day <- metadata_df$vaccination_day[1] # this will not work with the trial modelling when everyone has a different vccination day
  
  # This is for use with output of fit_rtss.R which is the infection records dataset 
  # and where all children are given the 3rd dose of the vaccine on day 0 of simulation 
  
  df2 <- df %>%
    filter(!is.na(detection_day)) %>%
    # filter(infectious_bite_day>0) %>% # to get rid of the build up 
    mutate(days_since_rtss = detection_day - vax_day,
           weeks_since_rtss = ceiling(days_since_rtss/7)) %>%
    group_by(arm, sim_id, weeks_since_rtss) %>%
    summarize(cases = n(),
              .groups = 'drop') %>%
    tidyr::complete(arm, sim_id, weeks_since_rtss = seq(min(weeks_since_rtss, na.rm = TRUE),
                                                        max(weeks_since_rtss, na.rm = TRUE), 1),
                    fill = list(cases = 0)) %>%
    mutate(pop = ifelse(arm == 'none', n_none, n_rtss)) %>% # I ran this with a pop of 1200 divided by rtss and no intervention 
    mutate(inci = cases / pop) %>% ungroup()
  
  # Filter so that it is at least 21 days since RTSS 
  df2 <- df2 %>%
    filter(weeks_since_rtss >= 3)
  
  # Extract no intervention  
  none <- df2 %>%
    filter(arm == 'none') %>% 
    rename(cases_none = cases, inci_none = inci, pop_none = pop) %>%
    select(-arm)
  
  # Extract RTSS cases 
  rtss <- df2 %>% 
    filter(arm == 'rtss')%>% 
    rename(cases_smc = cases, inci_rtss = inci, pop_rtss = pop) %>%
    select(-arm)
  
  d <- left_join(none, rtss, by = c('weeks_since_rtss','sim_id'))
  d$efficacy <-  1 - (d$inci_rtss / d$inci_none)
  
  return(d)
  
}

# # Function to calculate RTSS efficacy based on cumulative proportion infected
calc_rtss_efficacy_cumul <- function(df, params_row){
  n_rtss <- nrow(metadata_df[metadata_df$arm == 'rtss',])
  n_none <- nrow(metadata_df[metadata_df$arm == 'none',])
  vax_day <- metadata_df$vaccination_day[1] # this will not work with the trial modelling when everyone has a different vccination day
  
  # get first infection 
  df_ <- df %>%
    filter(!is.na(detection_day)) %>%
    arrange(detection_day, arm, rid) %>%
    group_by(arm, rid) %>%
    slice_min(detection_day) 
  
  df2 <- df_ %>%
    filter(!is.na(detection_day)) %>%
    # filter(infectious_bite_day>0) %>% # to get rid of the build up 
    mutate(days_since_rtss = detection_day - vax_day,
           weeks_since_rtss = ceiling(days_since_rtss/7)) %>%
    group_by(arm, sim_id, weeks_since_rtss) %>%
    summarize(cases = n(),
              .groups = 'drop') %>%
    tidyr::complete(arm, sim_id, weeks_since_rtss = seq(min(weeks_since_rtss, na.rm = TRUE),
                                                        max(weeks_since_rtss, na.rm = TRUE), 1),
                    fill = list(cases = 0)) %>%
    ungroup() %>% group_by(arm) %>%
    mutate(pop = ifelse(arm == 'none', n_none, n_rtss)) %>%
    mutate(cumulcases = cumsum(cases),
           cumulprop = cumulcases / pop) %>%
    ungroup() %>% select(-cases, -pop)
  
  # Extract no intervention  
  none <- df2 %>%
    filter(arm == 'none') %>% 
    rename(cumulprop_none = cumulprop,
           cumulcases_none = cumulcases) %>%
    select(-arm)
  
  # Extract SMC cases 
  rtss <- df2 %>% 
    filter(arm == 'rtss')%>% 
    rename(cumulprop_rtss = cumulprop,
           cumulcases_rtss = cumulcases) %>%
    select(-arm)
  
  d <- left_join(none, rtss, by = c('sim_id','weeks_since_rtss'))
  d$efficacy <-  1 - (d$cumulprop_rtss / d$cumulprop_none)
  
  return(d)
  
}


# # Function to calculate SMC efficacy based on Thompson et al. (cumulative proportion)
calc_smc_efficacy_cumul <- function(df, params_row, by_week = TRUE){
  smc_dose_days_ <- unlist(params_row$smc_dose_days)
  n_smc <- nrow(metadata_df[metadata_df$arm == 'smc',])
  n_none <- nrow(metadata_df[metadata_df$arm == 'none',])
  
  # get first infection 
  df_ <- df %>%
    filter(!is.na(detection_day)) %>%
    arrange(detection_day, arm, rid) %>%
    group_by(arm, rid) %>%
    slice_min(detection_day) 
  
  df_ <- df_ %>%
    mutate(most_recent_smc = sapply(detection_day, function(d){
      prior = smc_dose_days_[smc_dose_days_ < d]
      if(length(prior) == 0) NA_real_ else max(prior)
    })) %>%
    mutate(days_since_smc = detection_day - most_recent_smc,
           weeks_since_smc = ceiling(days_since_smc / 7)
    )
  
  if(by_week){
    dfweek <- df_ %>%
      group_by(arm, weeks_since_smc) %>%
      arrange(arm, weeks_since_smc) %>%
      summarise(cases = n(),
                .groups = 'drop') %>%
      tidyr::complete(arm, weeks_since_smc = seq(min(weeks_since_smc, na.rm = TRUE), 
                                                 max(weeks_since_smc, na.rm = TRUE), 1), 
                      fill = list(cases = 0)) %>%
      ungroup() %>% group_by(arm) %>%
      mutate(pop = ifelse(arm == 'none', n_none, n_smc)) %>%
      mutate(cumulcases = cumsum(cases),
             cumulprop = cumulcases / pop) %>%
      ungroup() 
    
    # Extract no intervention  
    none <- dfweek %>%
      filter(arm == 'none') %>% 
      rename(cumulprop_none = cumulprop,
             cumulcases_none = cumulcases) %>%
      select(-arm)
    
    # Extract SMC cases 
    smc <- dfweek %>% 
      filter(arm == 'smc')%>% 
      rename(cumulprop_smc = cumulprop,
             cumulcases_smc = cumulcases) %>%
      select(-arm)
    
    d <- left_join(none, smc, by = 'weeks_since_smc')
    d$efficacy <-  1 - (d$cumulprop_smc / d$cumulprop_none)
  }
  
  if(!by_week){
    dfday <- df_ %>%
      group_by(arm, days_since_smc) %>%
      arrange(arm, days_since_smc) %>%
      summarise(cases = n(),
                .groups = 'drop') %>%
      tidyr::complete(arm, days_since_smc = seq(min(days_since_smc, na.rm = TRUE), 
                                                max(days_since_smc, na.rm = TRUE), 1), 
                      fill = list(cases = 0)) %>%
      ungroup() %>% group_by(arm) %>%
      mutate(pop = ifelse(arm == 'none', n_none, n_smc)) %>%
      mutate(cumulcases = cumsum(cases),
             cumulprop = cumulcases / pop) %>%
      ungroup()
    
    # Extract no intervention  
    noneday <- dfday %>%
      filter(arm == 'none') %>% 
      rename(cumulprop_none = cumulprop,
             cumulcases_none = cumulcases) %>%
      select(-arm)
    
    # Extract SMC cases 
    smcday <- dfday %>% 
      filter(arm == 'smc')%>% 
      rename(cumulprop_smc = cumulprop,
             cumulcases_smc = cumulcases) %>%
      select(-arm)
    
    d <- left_join(noneday, smcday, 'days_since_smc')
    d$efficacy <-  1 - (d$cumulprop_smc / d$cumulprop_none)
  }
  
  return(d)
  
}

# function to calculate efficacy of SMC versus no intervention -- incidence based
calc_smc_efficacy <- function(df, params_row, by_week = TRUE){
  
  smc_dose_days_ <- unlist(params_row$smc_dose_days)
  n_smc <- nrow(metadata_df[metadata_df$arm == 'smc',])
  n_none <- nrow(metadata_df[metadata_df$arm == 'none',])
  
  if(by_week){
    df2 <- df %>%
      filter(!is.na(detection_day))  %>%# keep only cases 
      # group_by(arm) %>%
      # find most recent smc day 
      mutate(most_recent_smc = sapply(detection_day, function(d){
        prior = smc_dose_days_[smc_dose_days_ < d]
        if(length(prior) == 0) NA_real_ else max(prior)
      })) %>%
      mutate(days_since_smc = detection_day - most_recent_smc,
             weeks_since_smc = ceiling(days_since_smc / 7)
      ) %>%
      group_by(arm, weeks_since_smc) %>%
      summarise(cases = n(),
                .groups = 'drop') %>%
      tidyr::complete(arm, weeks_since_smc = seq(min(weeks_since_smc, na.rm = TRUE), 
                                                 max(weeks_since_smc, na.rm = TRUE), 1), 
                      fill = list(cases = 0)) %>%
      mutate(pop = ifelse(arm == 'none', n_none, n_smc)) %>%
      mutate(inci = cases / pop) %>% ungroup() %>%
      mutate(cases = replace_na(cases, 0),
             inci = replace_na(inci, 0))
    
    # Extract no intervention  
    none <- df2 %>%
      filter(arm == 'none') %>% 
      rename(cases_none = cases, inci_none = inci, pop_none = pop) %>%
      select(-arm)
    
    # Extract SMC cases 
    smc <- df2 %>% 
      filter(arm == 'smc')%>% 
      rename(cases_smc = cases, inci_smc = inci, pop_smc = pop) %>%
      select(-arm)
    
    d <- left_join(none, smc, by = 'weeks_since_smc')
    d$efficacy <-  1 - (d$inci_smc / d$inci_none)
  }
  
  if(!by_week){
    df2 <- df %>%
      filter(!is.na(detection_day))  %>%# keep only cases 
      # find most recent smc day 
      mutate(most_recent_smc = sapply(detection_day, function(d){
        prior = smc_dose_days_[smc_dose_days_ < d]
        if(length(prior) == 0) NA_real_ else max(prior)
      })) %>%
      mutate(days_since_smc = detection_day - most_recent_smc) %>%
      group_by(arm, days_since_smc) %>%
      summarise(cases = n(),
                .groups = 'drop') %>%
      tidyr::complete(arm, days_since_smc = seq(min(days_since_smc, na.rm = TRUE), 
                                                max(days_since_smc, na.rm = TRUE), 1), 
                      fill = list(cases = 0)) %>%
      mutate(pop = ifelse(arm == 'none', n_none, n_smc)) %>%
      mutate(inci = cases / pop) %>% ungroup() %>%
      mutate(cases = replace_na(cases, 0),
             inci = replace_na(inci, 0))
    
    # Extract no intervention  
    none <- df2 %>%
      filter(arm == 'none') %>% 
      rename(cases_none = cases, inci_none = inci, pop_none = pop) %>%
      select(-arm)
    
    # Extract SMC cases 
    smc <- df2 %>% 
      filter(arm == 'smc')%>% 
      rename(cases_smc = cases, inci_smc = inci, pop_smc = pop) %>%
      select(-arm)
    
    d <- left_join(none, smc, by = 'days_since_smc')
    d$efficacy <-  1 - (d$inci_smc / d$inci_none)
  }
  
  return(d)
}


