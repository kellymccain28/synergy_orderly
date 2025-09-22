# functions to do parameter sweeps 

# Main cohort simulation function
run_cohort_simulation <- function(metadata_df, 
                                  params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                                  base_inputs, 
                                  output_dir = 'simulation_outputs',
                                  allow_superinfections = TRUE, 
                                  save_outputs = TRUE) {
  
  # Extract parameters from the row
  max_SMC_kill_rate <- params_row$max_SMC_kill_rate
  smc_lambda <- params_row$lambda #17
  smc_kappa <- params_row$kappa #0.28 - these are default values (can't rmember where I got them though..)
  # SMC_decay <- params_row$SMC_decay
  # season_start_day <- params_row$season_start_day # always the same, but now with specific smc timings don't need this 
  # season_length <- 30*5#params_row$season_length
  lag_p_bite <- params_row$lag_p_bite
  # get specific probability of bite for each run 
  p_bite <- unlist(params_row$p_bite)
  
  message("Running simulation: ", params_row$sim_id)
  
  # Extract base inputs (these stay constant across parameter sweeps)
  trial_ts <- base_inputs$trial_timesteps#[[i]]
  burnin <- base_inputs$burnin
  threshold <- base_inputs$threshold
  VB <- base_inputs$VB
  tstep <- base_inputs$tstep
  t_liverstage <- base_inputs$t_liverstage
  weights <- base_inputs$weights#[[1]]
  # country <- base_inputs$country
  country_to_run <- base_inputs$country
  country_short <- str_sub(country_to_run, 1, 1)
  
  # now, just do 100 timesteps (200 days) for each infection, which is plenty of time 
  tt <- seq(0, 100, by = tstep)#seq(0, ((trial_ts + burnin) - t)/divide, by = tstep)
  
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
    country = character()
  )
  
  N <- nrow(metadata_df)
  
  # Make weights to prevent homogenoeous transmission 
  weights <- rexp(N) # Generate random weights to sample children at different probabilities (could also use rpois(N) or rexp(N))
  weights <- weights / sum(weights) # normalize to sum to 1
  
  message('calculating time since smc and parasite kill rate vector')
  # Calculate SMC kill rate vectors for each child  - not sure where these are from 
  # this is time since smc where the day is relative to the first day of the simulation (not including burnin)
  metadata_df$time_since_smc <- lapply(metadata_df$smc_dose_days, 
                                       calc_time_since_dose, 
                                       days = 0:(trial_ts+burnin))
  metadata_df$smckillvec <- lapply(metadata_df$time_since_smc,
                       function(.){
                       kill <- max_SMC_kill_rate * exp(-(./ smc_lambda)^smc_kappa)  # calculate kill rate with hill function
                       kill <- ifelse(is.na(kill),0, kill)                          # change NAs (when there is no SMC) to 0 
                       kill <- c(rep(0, burnin), kill)                              # pad front of vector with 0s for burnin period so that SMC begins after burnin
                       kill <- kill[seq_along(kill) %% 2 == 1]                      # keep only odd indices since timesteps in the within-host model are 2 days
                         })
  message('Calculated smc kill rate vectors for each child')
  # this kill vector has 25 2-day timesteps for burnin=50, then 548 2-day timesteps for the 1095 day simulation  
  
  
  # Initialize parasitemia storage
  parasitemia_storage <- vector("list", length = trial_ts + burnin)
  
  susceptibles <- rep(TRUE, N)
  
  # Choose probability of biting for baseline transmission/seasonality depending on country 
  # p_bite <- if(country_to_run == 'BF') p_bite_bfa else p_bite_mli
  
  # Main simulation loop 
  for (t in 1:(trial_ts + burnin)) {
    # message(str_glue("infection records has {nrow(infection_records)} entries"))
    
    # Determine number of new infectious bites on time t
    n_infectious_bites <- rpois(1, lambda = N * p_bite[t])
    
    if (n_infectious_bites == 0) next
    
    # Sample random children to be bitten 
    bites <- sample(1:N, size = n_infectious_bites, prob = weights, replace = TRUE)
    # message('There are ', length(bites), " mosquito bites on time ", t, ' out of ', (trial_ts+burnin))
    
    bit_kids <- unique(bites)
    # Remove any bit kids with no vaccine follow-up data 
    bit_kids <- bit_kids[bit_kids %in% metadata_df$rid] 
    message('There are ', length(bit_kids), " unique kids bitten on time ", t, ' out of ', (trial_ts+burnin))
    
    
    # Update susceptibility vector so that recovered kids would be susceptible again
    if (exists("infection_records") && nrow(infection_records) > 0 & !allow_superinfections) {
      
      # Filter to only records with actual detections
      detected_records <- infection_records[!is.na(infection_records$detection_day), ]
      
      if (nrow(detected_records) > 0) {
        # Get the most recent detection for each person
        most_recent <- detected_records %>%
          group_by(rid) %>%
          slice_max(detection_day, n = 1) %>%  # Get most recent detection_day for each rid
          ungroup()
        
        # Find who should recover today (day t)
        recovery_today <- (most_recent$recovery_day + burnin) == t # add burnin because the recovery day calculation is the external time but for this siulation we want ot know if on time t they could be infected
        # if the recovery day equals today then these are the ones that are newly susceptible 
        if (any(recovery_today)) {
          recovering_kids <- most_recent$rid[recovery_today]
          susceptibles[recovering_kids] <- TRUE
        }

      }
      }
    
    # Avoid superinfection -- if a kid is bitten and they aren't susceptible they are not added to the list of bit_kids
    bit_kids_t <- c()
    j <- 1
    for (kid in bit_kids){
      if (susceptibles[kid]){ 
        bit_kids_t[j] <- kid # Add this kid to the filtered list if the kid is susceptible
        
        j <- j + 1 # move to the next position in the filtered list 
      }
    }
    message(str_glue('{length(bit_kids_t)} / {length(bit_kids)} bit kids were susceptible'))
    bit_kids <- bit_kids_t
    
    # print(table(susceptibles))
    
    if(length(bit_kids) > 0) {
      # Estimate antibody level for each bitten child 
      if(t < burnin) {
        # No interventions during burnin
        PEV_vec <- rep(0, length(bit_kids))
        SMC_vec <- rep(0, length(bit_kids))
        # these both need a placeholder value (won't matter if PEV = 0 or SMC = 0)
        t_since_vax_vec <- rep(0, length(bit_kids))
        t_toboost1_vec <- rep(500, length(bit_kids))
        t_toboost2_vec <- rep(1000, length(bit_kids))
        
        SMC_kill_vec <- rep(list(rep(0, floor(trial_ts/2 + burnin))), length(bit_kids))
        SMC_timev <- rep(list(0:(length(SMC_kill_vec[[1]])-1)), length = length(bit_kids))
      } else {
        
        kid_metadata <- metadata_df[metadata_df$rid %in% bit_kids, ]
        
        # Create named vectors of PEVand vaccination days 
        pev_lookup <- setNames(kid_metadata$PEV, kid_metadata$rid)
        vax_day_lookup <- setNames(kid_metadata$vaccination_day, kid_metadata$rid)
        t_toboost1_lookup <- setNames(kid_metadata$t_to_boost1, kid_metadata$rid)
        t_toboost2_lookup <- setNames(kid_metadata$t_to_boost2, kid_metadata$rid)
        
        # Map to bit_kids order
        PEV_vec <- pev_lookup[as.character(bit_kids)]
        # Find the time since vaccination -- time between 3rd dose (vax day) and infectious bite
        t_since_vax_vec <- (t - burnin) - vax_day_lookup[as.character(bit_kids)] # here, a (-) value of vaxdaylookup indicates 
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
        # subset the kill rate vector to be from this external time 
        # divide max time by two because the vector is every two days 
        SMC_kill_vec <- lapply(SMC_kill_vec, function(smcvec){
          subset <- smcvec[floor((t + burnin) / 2) :length(smcvec)]
          return(subset)
        })
        
        SMC_timev <- rep(list(0:(length(SMC_kill_vec[[1]])-1)), length = length(bit_kids))#rep(list(0:(trial_ts/2 + burnin)), length = length(bit_kids))
        
        # Print intervention status
        message('PEV:', PEV_vec)
        message('SMC:', SMC_vec)
        }
      
      
      # Run within-host simulation with current parameters
      params_df <- data.frame(
        PEV_on = PEV_vec,
        SMC_on = SMC_vec,
        t_inf = t_since_vax_vec,
        infection_start_day = t,
        # season_start_day = season_start_day,  # Parameter from sweep
        # season_length = season_length,        # Parameter from sweep
        # smc_interval = smc_interval,          # Parameter from sweep
        SMC_time = I(SMC_timev),
        SMC_kill_vec = I(SMC_kill_vec), 
        rid = bit_kids,
        t_toboost1 = t_toboost1_vec,
        t_toboost2 = t_toboost2_vec
      )
      
      outputs <- pmap(params_df, 
                      function(PEV_on, 
                               SMC_on, 
                               t_inf, 
                               infection_start_day,
                               # season_start_day, 
                               # season_length, 
                               # smc_interval, 
                               SMC_time,
                               SMC_kill_vec,
                               rid,
                               t_toboost1,
                               t_toboost2) {
                        result <- run_process_model(
                          PEV_on = PEV_on,
                          SMC_on = SMC_on,
                          t_inf = t_inf,
                          VB = VB,
                          tt = tt,
                          max_SMC_kill_rate = max_SMC_kill_rate,  # Parameter from sweep
                          SMC_decay = SMC_decay,                  # Parameter from sweep
                          infection_start_day = infection_start_day,# external time that infection begins 
                          # season_start_day = season_start_day,# day of season start relative to Jan 1 of external year 
                          # season_length = season_length,# ~4 months (120 days)
                          # smc_interval = smc_interval,# how often are SMC rounds (days)
                          SMC_time = unlist(SMC_time),
                          SMC_kill_vec = unlist(SMC_kill_vec),
                          tboost1 = t_toboost1,
                          tboost2 = t_toboost2
                        )
                        
                        result$trajectory$rid <- rid
                        
                        return(result)
                      })
      
      # Vectorized data frame creation of infection records
      new_records <- data.frame(
        rid = bit_kids,
        time_ext = rep(t - burnin, length(bit_kids)),                                                      # external cohort time (t is the cohort time, then we wnat to scale to be + if after burnin)
        infectious_bite_day = rep((t - burnin) - t_liverstage, length(bit_kids)),                          # bitten on day t, then assuming the liver stage takes t_liverstage days
        BSinfection_day = rep((t - burnin), length(bit_kids)),                                             # after liver stage, the BS begins #+ t_liverstage
        threshold_day = sapply(outputs, function(x) x$threshold_day),               # days since BS starts that threshold is reached
        detection_day = (t - burnin) + sapply(outputs, function(x) x$threshold_day),# day in cohort simulation that threshold is reached, threshold is the day since BS infection that reaches threshold #+ t_liverstage 
        t_toreach_threshold = sapply(outputs, function(x) x$threshold_day) + t_liverstage,   # time to reach threshold value / detection since the bite  
        vaccination_day = if(t < burnin) rep(NA, length(bit_kids)) else kid_metadata$vaccination_day,   # day of vaccination relative to the start of follow-up (day 0 external time)
        prob_bite = rep(p_bite[t], length(bit_kids)),
        recovery_day = ((t - burnin) + sapply(outputs, function(x) x$threshold_day)) + 12 - t_liverstage, # day that the child would be 'recovered' if we assume that a child is treated and has a period of prophylaxis for 12 days after detection day 
                                                                      # (90% at 12 days in paper but here, assuming 100% for 12 days) after the day of treatment and that all infectiosn are treated with AL 10.1038/ncomms6606
        country = country_to_run
      ) 
      
      infection_records <- rbind(infection_records, new_records)
      
      # Only update the susceptibility vector if we are not allowing superinfections 
      if(!allow_superinfections){
        # mark as non-susceptible only those who have a valid detection_day
        detected_kids <- bit_kids[!is.na(new_records$detection_day)]
        # print(detected_kids)
        if (length(detected_kids) > 0) {
          susceptibles[detected_kids] <- FALSE
        } 
      }
      
      # Vectorized parasitemia storage creation
      parasitemia_data <- map2(outputs, bit_kids, function(output, kid) {
        
        output$trajectory %>%
          mutate(
            day1_BSinfection = t - burnin,#+ t_liverstage 
            detection_day = t + (output$threshold_day) - burnin,# do not need to multiply threshold day by 2 here since already done in run_process_model #+ t_liverstage 
            time_ext = time*2 + (t - burnin),# + infection_start_day/2,# external time should be dependent on when infection was, relative to external time(inf_start_day aka t); time (model time) has already been multiplied by 2 and related to infection time  #time*2 + (t - 1) - burnin,#+ t_liverstage
            arm = metadata_df[metadata_df$rid == kid, ]$arm,
            t = t
          )
      })
      
      # Store all parasitemia data for this time step
      parasitemia_storage[[t]] <- parasitemia_data
    }
  }
  
  message("finished sim")
  
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
  
  # Process infection records
  child_counts <- metadata_df %>%
    distinct(rid, arm, country) %>%
    count(arm, name = "children_in_group")
  
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
  
  # Create diagnostic plots
  plots <- create_diagnostic_plots(parasitemia_df, infection_records_final,
                                   max_SMC_kill_rate,
                                   smc_lambda, smc_kappa)
  message('made plots')
  
  # Prepare results
  results <- list(
    sim_id = params_row$sim_id,
    parameters = params_row,
    infection_records = infection_records_final,
    child_metadata = metadata_df,
    parasitemia_data = parasitemia_df,
    diagnostic_plots = plots
  )
  
  message('prepared results')

    # Save outputs if requested
  if (save_outputs && !is.null(output_dir)) {
    
    save_simulation_outputs(results, output_dir)
    message('saved outputs')
  }
  
  return(results)
}


# Function to create diagnostic plots
create_diagnostic_plots <- function(parasitemia_df, 
                                    infection_records_final,
                                    max_SMC_kill_rate, smc_lambda, smc_kappa) {
  
  # Plot 1: Parasitemia over time
  p1 <- ggplot(parasitemia_df) +
    geom_line(aes(x=time_ext/365.25, y = parasites, 
                  group = as.factor(paste0(rid, ", day ", day1_BSinfection)), 
                  color = det), alpha = 0.5) +
    scale_y_log10() +
    geom_hline(aes(yintercept = 10), linetype = 2, color = 'darkred', linewidth = 1) +
    geom_hline(aes(yintercept = 1e-5), linetype = 2, color = 'darkgreen', linewidth = 1) +
    geom_vline(aes(xintercept = 0), color = 'darkturquoise', linetype = 2, linewidth = 1) +
    facet_wrap(~arm) +
    theme_bw() +
    theme(legend.position = 'none') +
    labs(title = "Parasitemia trajectories",
         caption = str_glue("SMC max kill rate: {max_SMC_kill_rate}, SMC lambda: {smc_lambda}, SMC kappa: {smc_kappa}"))
  
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
    parasitemia_trajectories = p1,
    proportion_detectable = p2,
    incidence = p3,
    cumulative_infections = p4,
    prob_infectious_bite = p5
  ))
}

# Function to save simulation outputs
save_simulation_outputs <- function(results, output_dir) {
  
  # dir.create('outputs/', showWarnings = FALSE)
  
  sim_dir <- file.path(output_dir, results$sim_id)
  
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir, recursive = TRUE)
  }
  
  plots_dir <- file.path(sim_dir, "/plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir)
  }
  
  # Save input files 
  # write.csv(params_df, file.path(sim_dir, "parameter_grid.csv"), row.names = FALSE)
  saveRDS(base_inputs, file.path(sim_dir, "base_inputs.rds"))
  
  # Save data files
  saveRDS(results$infection_records, file.path(sim_dir, "infection_records.rds"))
  saveRDS(results$child_metadata, file.path(sim_dir, "child_metadata.rds"))
  saveRDS(results$parasitemia_data, file.path(sim_dir, "parasitemia_data.rds"))
  saveRDS(results$parameters, file.path(sim_dir, "parameters.rds"))
  
  # Save plots
  ggsave(file.path(plots_dir, "parasitemia_trajectories.png"), 
         results$diagnostic_plots$parasitemia_trajectories, 
         height = 8, width = 12)
  ggsave(file.path(plots_dir, "proportion_detectable.png"), 
         results$diagnostic_plots$proportion_detectable, 
         height = 8, width = 12)
  ggsave(file.path(plots_dir, "incidence.png"), 
         results$diagnostic_plots$incidence, 
         height = 6, width = 14)
  ggsave(file.path(plots_dir, "cumulative_infections.png"), 
         results$diagnostic_plots$cumulative_infections, 
         height = 8, width = 12)
  ggsave(file.path(plots_dir, "prob_infectious_bite.png"), 
         results$diagnostic_plots$prob_infectious_bite, 
         height = 8, width = 12)
  
  message("Saved outputs for ", results$sim_id, " to ", sim_dir)
}

# Main function to run parameter sweep
run_parameter_sweep <- function(metadata_df, 
                                parameters_df, 
                                base_inputs, 
                                output_dir = "simulation_outputs", 
                                save_outputs = TRUE,
                                parallel = FALSE, n_cores = NULL,
                                allow_superinfections = TRUE) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  parameters_df$sim_id <- paste0('parameter_set_',rownames(parameters_df),"_", country_to_run, "_",allow_superinfections)#paste0("sim_", "SMCkill", parameters_df$max_SMC_kill_rate,"_SMCdecay", parameters_df$SMC_decay,"_season", parameters_df$season_start_day)
  
  # Save parameter grid
  saveRDS(parameters_df, file.path(output_dir, "parameter_grid.rds"))
  saveRDS(base_inputs, file.path(output_dir, "base_inputs.rds"))
  cyphr::encrypt(saveRDS(metadata_df, file.path(output_dir, "metadata_df.rds")), key)
  
  if (parallel) {
    # if (is.null(n_cores)) {
    #   n_cores <- parallel::detectCores() - 1
    # }
    
    message('Running in parallel is not working at the moment')
    # message("Running ", nrow(params_df), " simulations in parallel using ", n_cores, " cores")
    
    # Set up parallel backend (works on Windows, Mac, Linux)
    # plan(multisession, workers = n_cores)
    
    # # Split params_df into list of rows
    # param_list <- split(params_df, seq(nrow(params_df)))
    # 
    # results <- future_lapply(param_list, 
    #                          function(param_row) {
    #                            run_cohort_simulation(param_row, base_inputs, output_dir,
    #                                                  save_outputs)
    #                          },
    #                          future.seed = TRUE)
    # Split parameter grid into batches
    # n_per_batch <- 3
    # param_batches <- split(params_df, ceiling(seq_len(nrow(params_df))/n_per_batch))
    # 
    # # Run each batch
    # results <- future_lapply(param_batches, 
    #                          function(param_row) {
    #   message("Running batch ", i, " of ", length(param_batches))
    #   run_parameter_sweep(param_batches[[i]], base_inputs)
    # })
    
    
    # Reset to sequential processing
    # plan(sequential)
  } else {
    message("Running ", nrow(params_df), " simulations sequentially")
    
    results2 <- list()
    for (i in 1:nrow(params_df)) {
      results2[[i]] <- run_cohort_simulation(metadata_df, parameters_df[i, ], base_inputs, output_dir,
                                            save_outputs = TRUE,
                                            allow_superinfections = allow_superinfections)#results[[i]] <- 
      
      message("done simulation ", i, " of ", nrow(parameters_df))
      
    }
  }
  
  return(list(
    results2
    # summary_results = summary_results
  ))
  
}
