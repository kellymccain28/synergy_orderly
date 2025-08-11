# functions to do parameter sweeps 

# Main cohort simulation function
run_cohort_simulation <- function(params_row, 
                                  base_inputs, 
                                  output_dir = 'simulation_outputs',
                                  save_outputs = TRUE) {
  
  # Extract parameters from the row
  max_SMC_kill_rate <- params_row$max_SMC_kill_rate
  SMC_decay <- params_row$SMC_decay
  season_start_day <- params_row$season_start_day
  season_length <- 30*5#params_row$season_length
  
  # Create unique simulation ID
  sim_id <- paste0("sim_", 
                   "SMCkill", max_SMC_kill_rate,
                   "_SMCdecay", SMC_decay,
                   "_season", season_start_day,
                   "_length", season_length)
  
  message("Running simulation: ", sim_id)
  
  # Extract base inputs (these stay constant across parameter sweeps)
  N <- base_inputs$N_children
  trial_ts <- base_inputs$trial_timesteps#[[i]]
  burnin <- base_inputs$burnin
  p_bite_bfa <- base_inputs$p_bite_bfa
  p_bite_mli <- base_inputs$p_bite_mli
  threshold <- base_inputs$threshold
  VB <- base_inputs$VB
  tstep <- base_inputs$tstep
  t_liverstage <- base_inputs$t_liverstage
  weights <- base_inputs$weights#[[1]]
  vax_day <- base_inputs$vax_day
  country <- base_inputs$country
  smc_interval <- base_inputs$smc_interval
  
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
    p_bite = double()
  )
  
  # Create child metadata 
  # when i simulate the real cohort, this will be more specific to the information in the data from the trial 
  metadata_child <- data.frame(
    rid = 1:N,
    vaccination_day = vax_day,
    PEV = rbinom(N, 1, 0.5),
    SMC = c(rep(1, N/2), rep(0, N/2))
  ) %>%
    mutate(arm = case_when(
      PEV == 1 & SMC == 1 ~ 'vaxsmc',
      PEV == 1 & SMC == 0 ~ 'vax',
      PEV == 0 & SMC == 1 ~ 'smc',
      TRUE ~ 'none'))
  
  # Initialize parasitemia storage
  parasitemia_storage <- vector("list", length = trial_ts + burnin)
  
  # Main simulation loop (adapted from your original)
  for (t in 1:(trial_ts + burnin)) {
    
    # Choose probability of biting for baseline transmission/seasonality depending on country 
    p_bite <- if(country == 'BF') p_bite_bfa else p_bite_mli
    
    # Determine number of new infectious bites on time t
    n_infectious_bites <- rpois(1, lambda = N * p_bite[t])
    
    if (n_infectious_bites == 0) next
    
    # Sample random children to be bitten 
    bites <- sample(1:N, size = n_infectious_bites, prob = weights, replace = TRUE)
    message('There are ', length(bites), " mosquito bites on time ", t, ' out of ', (trial_ts+burnin))
    
    bit_kids <- unique(bites)
    message('There are ', length(bit_kids), " unique kids bitten on time ", t, ' out of ', (trial_ts+burnin))
    
    if(length(bit_kids) > 0) {
      # Estimate antibody level for each bitten child 
      if(t < burnin) {
        # No interventions during burnin
        PEV_vec <- rep(0, length(bit_kids))
        SMC_vec <- rep(0, length(bit_kids))
        # these both need a placeholder value (won't matter if PEV = 0 or SMC = 0)
        t_since_vax_vec <- rep(0, length(bit_kids))
      } else {
        
        # Vectorized intervention calculations
        kid_metadata <- metadata_child[metadata_child$rid %in% bit_kids, ]
        
        # Create named vectors
        pev_lookup <- setNames(kid_metadata$PEV, kid_metadata$rid)
        smc_lookup <- setNames(kid_metadata$SMC, kid_metadata$rid)
        vax_day_lookup <- setNames(kid_metadata$vaccination_day, kid_metadata$rid)
        
        # Map to bit_kids order
        PEV_vec <- pev_lookup[as.character(bit_kids)]
        SMC_vec <- smc_lookup[as.character(bit_kids)]
        t_since_vax_vec <- (t - burnin) - vax_day_lookup[as.character(bit_kids)]
        
        # Get intervention status
        message('PEV:', PEV_vec)
        message('SMC:', SMC_vec)
        }
      
      # Run within-host simulation with current parameters
      params_df <- data.frame(
        PEV_on = PEV_vec,
        SMC_on = SMC_vec,
        t_inf = t_since_vax_vec,
        infection_start_day = t,
        season_start_day = season_start_day,  # Parameter from sweep
        season_length = season_length,        # Parameter from sweep
        smc_interval = smc_interval,          # Parameter from sweep
        rid = bit_kids
      )
      
      # Find time length until end of follow-up to run infection sim  
      tt_until_end_cohort <- seq(0, ((trial_ts + burnin) - t)/divide, by = tstep)
      
      outputs <- pmap(params_df, 
                      function(PEV_on, 
                               SMC_on, 
                               t_inf, 
                               infection_start_day,
                               season_start_day, 
                               season_length, 
                               smc_interval, 
                               rid) {
                        result <- run_process_model(
                          PEV_on = PEV_on,
                          SMC_on = SMC_on,
                          t_inf = t_inf,
                          VB = VB,
                          tt = tt_until_end_cohort,
                          max_SMC_kill_rate = max_SMC_kill_rate,  # Parameter from sweep
                          SMC_decay = SMC_decay,                  # Parameter from sweep
                          infection_start_day = infection_start_day,# external time that infection begins 
                          season_start_day = season_start_day,# day of season start relative to Jan 1 of external year 
                          season_length = season_length,# ~4 months (120 days)
                          smc_interval = smc_interval# how often are SMC rounds (days)
                        )
                        
                        result$trajectory$rid <- rid
                        
                        return(result)
                      })
      
      # Vectorized data frame creation for infection records
      new_records <- data.frame(
        time_ext = t - burnin,                                                      # external cohort time (t is the cohort time, then we wnat to scale to be + if after burnin)
        rid = bit_kids,                                                        # number of bit kids
        infectious_bite_day = (t - burnin) - t_liverstage,                          # bitten on day t, then assuming the liver stage takes t_liverstage days
        BSinfection_day = (t - burnin),                                             # after liver stage, the BS begins #+ t_liverstage
        threshold_day = sapply(outputs, function(x) x$threshold_day),               # days since BS starts that threshold is reached
        detection_day = (t - burnin) + sapply(outputs, function(x) x$threshold_day),# day in cohort simulation that threshold is reached, threshold is the day since BS infection that reaches threshold #+ t_liverstage 
        t_toreach_threshold = sapply(outputs, function(x) x$threshold_day) + t_liverstage,   # time to reach threshold value / detection since the bite  
        vaccination_day = if(t < burnin) rep(vax_day, length(bit_kids)) else kid_metadata$vaccination_day,   # day of vaccination relative to the start of follow-up (day 0 external time)
        prob_bite = p_bite[t]
      ) 
      
      infection_records <- rbind(infection_records, new_records)
      
      # Vectorized parasitemia storage creation
      parasitemia_data <- map2(outputs, bit_kids, function(output, kid) {
        
        output$trajectory %>%
          mutate(
            day1_BSinfection = t - burnin,#+ t_liverstage 
            detection_day = t + (output$threshold_day) - burnin,# do not need to multiply threshold day by 2 here since already done in run_process_model #+ t_liverstage 
            time_ext = time*2 + (t - burnin),# + infection_start_day/2,# external time should be dependent on when infection was, relative to external time(inf_start_day aka t); time (model time) has already been multiplied by 2 and related to infection time  #time*2 + (t - 1) - burnin,#+ t_liverstage
            arm = metadata_child[metadata_child$rid == kid, ]$arm,
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
    filter(detection_day > 0 | is.na(detection_day)) %>%
    group_by(rid, day1_BSinfection) %>%
    mutate(
    #   child_dayinf = paste0(rid, ", day ", day1_BSinfection),
      det = ifelse(!is.na(detection_day), 1, 0)
    ) 
    
    # keep a small sample of the parasitemia dataset
    parasitemia_df <- parasitemia_df[parasitemia_df$rid %in% 
                     (parasitemia_df %>% 
                        distinct(rid, arm) %>% 
                        group_by(arm) %>% 
                        slice_sample(n=50) %>% 
                        pull(rid)), ]
  
  # Process infection records
  child_counts <- metadata_child %>%
    distinct(rid, arm) %>%
    count(arm, name = "children_in_group")
  
  infection_records_final <- infection_records %>%
    filter(detection_day > 0 | is.na(detection_day)) %>%
    left_join(metadata_child %>% select(-vaccination_day), by = 'rid') %>%
    left_join(child_counts, by = 'arm') %>%
    mutate(detectable = ifelse(is.na(threshold_day), 0, 1)) %>%
    group_by(arm) %>%
    arrange(detection_day) %>%
    mutate(cumul_inf = cumsum(detectable))
  
  message('processed sim')
  
  # Create diagnostic plots
  plots <- create_diagnostic_plots(parasitemia_df, infection_records_final,
                                   max_SMC_kill_rate, SMC_decay, season_start_day, season_length)
  message('made plots')
  
  # Prepare results
  results <- list(
    sim_id = sim_id,
    parameters = as.list(params_row),
    infection_records = infection_records_final,
    child_metadata = metadata_child,
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
                                    max_SMC_kill_rate, SMC_decay, season_start_day, season_length) {
  
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
         caption = str_glue("SMC max kill rate: {max_SMC_kill_rate}, SMC decay: {SMC_decay}\nseason start: {season_start_day}, season length: {season_length}"))
  
  # Plot 2: Proportion detectable
  p2 <- ggplot(infection_records_final) +
    geom_bar(aes(x = arm, group = as.factor(detectable), 
                 fill = as.factor(detectable)), position = 'fill') +
    scale_fill_manual(values = c('darkmagenta','goldenrod')) +
    labs(fill = 'Detectable', title = "Proportion detectable infections",
         caption = str_glue("SMC max kill rate: {max_SMC_kill_rate}, SMC decay: {SMC_decay}, season start: {season_start_day},\nseason length: {season_length}"))
  
  # Plot 3: Incidence
  p3 <- ggplot(infection_records_final) +
    geom_line(aes(x = detection_day, color=arm), stat = 'count') +
    facet_wrap(~arm) +
    labs(y = 'N infections', x = 'Day since start of follow up period',
         title = "Incidence over time") +
    theme(legend.position = 'none',
          caption = str_glue("SMC max kill rate: {max_SMC_kill_rate}, SMC decay: {SMC_decay}, season start: {season_start_day},\nseason length: {season_length}"))
  
  # Plot 4: Cumulative infections
  p4 <- ggplot(infection_records_final) +
    geom_line(aes(x = detection_day, y = cumul_inf, color=arm)) +
    labs(y = 'N infections', 
         x = 'Day since start of follow up period',
         title = "Cumulative infections",
         caption = str_glue("SMC max kill rate: {max_SMC_kill_rate}, SMC decay: {SMC_decay}, season start: {season_start_day},\nseason length: {season_length}"))
  
  return(list(
    parasitemia_trajectories = p1,
    proportion_detectable = p2,
    incidence = p3,
    cumulative_infections = p4
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
  # saveRDS(base_inputs, file.path(sim_dir, "base_inputs.rds"))
  
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
         height = 8, width = 12)
  ggsave(file.path(plots_dir, "cumulative_infections.png"), 
         results$diagnostic_plots$cumulative_infections, 
         height = 8, width = 12)
  
  message("Saved outputs for ", results$sim_id, " to ", sim_dir)
}

# Main function to run parameter sweep
run_parameter_sweep <- function(params_df, base_inputs, output_dir = "simulation_outputs", 
                                save_outputs = TRUE,
                                parallel = FALSE, n_cores = NULL) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save parameter grid
  write.csv(params_df, file.path(output_dir, "parameter_grid.csv"), row.names = FALSE)
  saveRDS(base_inputs, file.path(output_dir, "base_inputs.rds"))
  
  
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
      results2[[i]] <- run_cohort_simulation(params_df[i, ], base_inputs, output_dir,
                                            save_outputs = TRUE)#results[[i]] <- 
      
      message("done simulation ", i, " of ", nrow(params_df))
      
      # Calculate log-likelihood and append it to with params_df
    }
  }
  
  # # Save summary results
  # summary_results <- map_dfr(results, function(x) {
  #   data.frame(
  #     sim_id = x$sim_id,
  #     x$parameters,
  #     n_detectable = sum(x$infection_records$detectable, na.rm = TRUE),
  #     prop_detectable = mean(x$infection_records$detectable, na.rm = TRUE)
  #   )
  # })
  # 
  # write.csv(summary_results, file.path(output_dir, "summary_results.csv"), row.names = FALSE)

  message("Parameter sweep complete! Results saved to: ", output_dir)
  
  return(list(
    results2
    # summary_results = summary_results
  ))
  
}
