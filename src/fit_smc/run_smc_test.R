# script to quickly test the smc parameters from optimization after running fit_smc.R

run_smc_test <- function(path = "R:/Kelly/synergy_orderly",
                        n_param_sets,
                        treatment_prob = 0.9, # default is 1 (which gives children prophylaxis)
                        N = 800){
  # Script to fit smc parameters to Hayley's curve 
  library(lhs)
  library(odin2)
  library(dust2)
  library(purrr)
  library(dplyr)
  library(mgcv)
  library(umbrella)
  library(lhs)
  library(cyphr)
  library(survival)
  library(broom)
  library(survminer)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  
  # Source antibody function
  source(paste0(path, "/shared/rtss.R"))
  # Source helper functions
  source(paste0(path, "/shared/helper_functions.R"))
  # Load the within-host model
  gen_bs <- odin2::odin(paste0(path, "/shared/smc_rtss.R"))
  # Source the utils functions
  source(paste0(path, "/shared/cohort_sim_utils.R"))
  # Source processing functions
  source(paste0(path, "/src/fit_smc/calculate_efficacy_likelihood.R"))
  
  trial_ts = 80# trial timesteps in cohort simulation
  # sim_allow_superinfections = TRUE # TRUE or FALSE
  country_to_run = 'generic'# BF or Mali, or if generic, then 'generic' which means that metadata_df is different.
  country_short = 'g'
  n_param_sets = n_param_sets
  N = N
  vax_day = 0 # unlike the model sim, this is in days (not timesteps)
  
  n_particles = 1L
  n_threads = 1L
  burnints = 30
  threshold = 5000
  tstep = 1
  t_liverstage = 8
  VB = 1e6
  divide = if(tstep == 1) 2 else 1
  
  treatment_probability = treatment_prob # in trial, everyone who was diagnosed with clincial malaria was treated 
  successful_treatment_probability = 0.9 # AL treatment protects up to 90% for 12 days SI of Commun. 5:5606 doi: 10.1038/ncomms6606 (this is what is in malsim)
  
  # Set up base inputs (these don't vary across parameter sweep)
  base_inputs <- list(
    trial_timesteps = trial_ts,
    burnin = burnints,
    threshold = threshold,
    VB = VB,
    tstep = tstep,
    t_liverstage = t_liverstage,
    country = country_to_run,
    country_short = 'g',
    treatment_probability = treatment_probability, 
    successful_treatment_probability = successful_treatment_probability
  )
  
  # Set up grid of parameters
  # param_ranges <- list(
  #   max_SMC_kill_rate = c(1, 3),# parasites per uL per 2-day timestep
  #   lambda = c(15, 19),
  #   kappa = c(0.1, 9)
  # )
  # # Generate LHS samples
  # A <- randomLHS(n_param_sets, 3)
  # # Scale to parameter ranges
  # params_df <- data.frame(
  #   max_SMC_kill_rate = qunif(A[,1], param_ranges$max_SMC_kill_rate[1], param_ranges$max_SMC_kill_rate[2]),
  #   lambda = qunif(A[,2], param_ranges$lambda[1], param_ranges$lambda[2]),
  #   kappa = qunif(A[,3], param_ranges$kappa[1], param_ranges$kappa[2])
  # )
  # "fitted" parameter values for SMC
  # first are the fitted parameters from 29 nov 
  # params_df <- data.frame(
  #   max_SMC_kill_rate = rep(c(2.333333, 2.333333, 2.333333, 1.888889, 2.333333, 
  #                         # now the best one from original optimization on 28th
  #                         2.447318,
  #                         # from old grid search 
  #                         2.071429), 5),
  #   lambda = rep(c(16.6667, 16.66667, 15, 25, 18.33333,
  #              # now the best one from original optimization on 28th
  #              17.06041, 
  #              # from old grid search 
  #              20.71429),5),
  #   kappa = rep(c(0.3444444, 0.2222222, 0.2222222, 0.5888889, 0.3444444,
  #             # now the best one from original optimization on 28th
  #             0.893555,
  #             # from old grid search 
  #             0.500000 ),5),
  #   sim_id = seq(1:35),
  #   repnum = rep(seq(1:7),5))
  # params_df <- data.frame(
  #   max_SMC_kill_rate = rep(c(2.333333, 2.333333, 2.333333), 10),
  #   lambda = rep(c(16.6667, 16.66667, 18.33333),10),
  #   kappa = rep(c(0.3444444, 0.2222222, 0.3444444),10),
  #   sim_id = seq(1:30),
  #   repnum = rep(seq(1:3),10))
  # params_df <- readRDS(paste0(path, '/src/fit_smc/outputs/top_20params_1201.rds'))
  # params_df$repnum <- rownames(params_df)
  # params_df <- rbind(params_df, params_df, params_df)
  # params_df$sim_id <- rownames(params_df)
  
  # Do multiple repetitionso f the final parameters 
  params_df <- data.frame(
    max_SMC_kill_rate = rep(2.333333, n_param_sets),
      lambda = rep(16.6667,n_param_sets),
      kappa = rep(0.2222222, n_param_sets),
      sim_id = seq(1:n_param_sets)) 
  
  # params_df$sim_id <- paste0('parameter_set_', rownames(params_df),"_", country_to_run, "_", treatment_probability)
  prob_bite_generic <- readRDS(paste0(path, '/archive/fit_rainfall/20251009-144330-1d355186/prob_bite_generic.rds'))
  prob_bite_generic$prob_infectious_bite = 0.3
  p_bitevector <- calc_lagged_vectors(prob_bite_generic, 0, burnints = burnints) # no lagged values
  
  params_df$lag_p_bite <- 0
  params_df$season_start_day <- 0
  
  # SMC delivery
  params_df <- params_df %>%
    rowwise() %>%
    mutate(smc_dose_days = 10
    ) %>%
    ungroup()
  
  parameters_df <- params_df %>%
    mutate(
      p_bite = purrr::map(lag_p_bite, ~p_bitevector[[paste0("lag_", .x)]])
    )
  
  # Make list of parameters instead of df
  params_list <- split(parameters_df, seq(nrow(parameters_df))) #%>%
  saveRDS(parameters_df, file.path(path, 'parameters_df.rds'))
  # Make metadata
  # Day of intervention (0 = start of follow-up; - values are before follow-up; + values after follow-up) - where 0 is also end of burnin
  # these get converted later to the correct directon - i.e. vaccine before follow-up will be +, smc before follow up will be -
  # vax_day is the 3rd primary dose (when we assume that efficacy begins)
  
  metadata_df <- data.frame(
    rid = 1:N,
    vaccination_day = vax_day,#sample(-1:0, N, replace = TRUE), the model takes abs(days before start of follow-up)
    PEV = 0,#rep(0,N),#
    SMC = c(rep(1, N/2), rep(0, N/2))#rbinom(N, 1, 0.5)
  ) %>%
    mutate(arm = case_when(
      PEV == 1 & SMC == 1 ~ 'both',
      PEV == 1 & SMC == 0 ~ 'rtss',
      PEV == 0 & SMC == 1 ~ 'smc',
      TRUE ~ 'none')) %>%
    mutate(t_to_boost1 = trial_ts ,
           t_to_boost2 = trial_ts + burnints,
           country = country_to_run) %>%
    mutate(rid_original = paste0(country_short, sprintf("%04d", rid)),
           country = 'generic',
           v1_date = as.Date('2017-04-01'))
  
  # saveRDS(metadata_df, "R:/Kelly/synergy_orderly/src/fit_smc/outputs/metadata_df.rds")
  
  
  # Make grid search list 
  # Define ranges
  # kill_vals  <- seq(2, 3,  length.out = 10)
  # lambda_vals <- seq(15, 19, length.out = 10)
  # kappa_vals  <- seq(0.1, 0.5,  length.out = 10)
  # 
  # # Cartesian product
  # grid <- expand.grid(
  #   max_SMC_kill_rate = kill_vals,
  #   lambda            = lambda_vals,
  #   kappa             = kappa_vals
  # ) %>%
  #   mutate(lag_p_bite = 0,
  #          smc_dose_days = 10,
  #          p_bite = parameters_df$p_bite[1])
  # params_list <- split(grid, seq(nrow(grid)))
  
  # Get the observed efficacy to compare to 
  observed_efficacy <- read.csv(paste0(path, '/shared/smc_fits_hayley.csv')) %>%
    mutate(weeks_since_smc = ceiling(day_since_smc / 7)) %>%
    group_by(weeks_since_smc) %>%
    mutate(efficacy = mean(efficacy))
  
  # Run simulation
  cluster_cores <- Sys.getenv("CCP_NUMCPUS")
  # if (cluster_cores == "") {
  #   cluster_cores <- 8
  # }
  
  if (cluster_cores == "") {
    message("running in serial (on a laptop?)")
    message("Running ", nrow(params_df), " simulations sequentially")
    
    # Grid search 
    grid_search_outputs <- lapply(
      params_list, 
      function(params_row){
        startsim <- Sys.time()
        o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                                   metadata_df,
                                   base_inputs,
                                   output_dir = 'R:/Kelly/src/fit_smc/simulation_outputs',
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
        
        # Remove NAs before calculating likelihood
        matched_complete <- matched %>%
          filter(!is.na(observed_efficacy), !is.na(predicted_efficacy))
        
        # Calculate mean least squares
        mls <- mean((matched_complete$observed_efficacy - matched_complete$predicted_efficacy)^2)
        
        return(list(mls = mls, 
                    efficacy = eff,
                    params = params_row))
      }
    )
    # 
    # results2 <- lapply(params_list,
    #                    function(params_row){
    #                      o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
    #                                                 metadata_df,
    #                                                 base_inputs,
    #                                                 output_dir = 'R:/Kelly/synergy_orderly/src/fit_smc/simulation_outputs/',
    #                                                 # allow_superinfections = TRUE,
    #                                                 return_parasitemia = TRUE,
    #                                                 save_outputs = FALSE)
    #                      message('finished simulation')
    #                      eff <- calc_smc_efficacy(o$infection_records,
    #                                               params_row,
    #                                               by_week = TRUE)
    #                      eff_daily <- calc_smc_efficacy(o$infection_records,
    #                                                     params_row,
    #                                                     by_week = FALSE)
    #                      eff$sim_id <- params_row$sim_id
    #                      eff_daily$sim_id <- params_row$sim_id
    #                      
    #                      # Efficacy with cumulative proportion
    #                      eff_cumul <- calc_smc_efficacy_cumul(o$infection_records,
    #                                                           params_row,
    #                                                           by_week = TRUE)
    #                      eff_daily_cumul <- calc_smc_efficacy_cumul(o$infection_records,
    #                                                                 params_row,
    #                                                                 by_week = FALSE)
    #                      eff_cumul$sim_id <- params_row$sim_id
    #                      eff_daily_cumul$sim_id <- params_row$sim_id
    #                      
    #                      return(list(efficacy_weekly = eff,
    #                                  efficacy_daily = eff_daily,
    #                                  efficacy_weekly_cumul = eff_cumul, 
    #                                  efficacy_daily_cumul = eff_daily_cumul,
    #                                  parasitemia = o$parasitemia_data %>% mutate(sim_id = params_row$sim_id),
    #                                  infection_records = o$infection_records %>% mutate(sim_id = params_row$sim_id),
    #                                  params = params_row))
    #                    })
    
  } else {
    message(sprintf("running in parallel on %s (on the cluster?)", cluster_cores))
    cl <- parallel::makeCluster(as.integer(cluster_cores),
                                outfile ="")
    invisible(parallel::clusterCall(cl, ".libPaths", .libPaths()))
    parallel::clusterCall(cl, function() {
      message('running')
      library(odin2)
      library(ggplot2)
      library(dust2)
      library(tidyverse)
      library(mgcv)
      library(umbrella)
      library(lhs)
      # library(orderly)
      library(retry)
      library(cyphr)
      library(survival)
      library(broom)
      library(survminer)
      library(tidyr)
      library(purrr)
      library(stringr)
      
      source('R:/Kelly/synergy_orderly/shared/cohort_sim_utils.R')
      source('R:/Kelly/synergy_orderly/shared/helper_functions.R')
      source("R:/Kelly/synergy_orderly/shared/rtss.R")
      source("R:/Kelly/synergy_orderly/shared/likelihood.R")
      source("R:/Kelly/synergy_orderly/src/fit_smc/calculate_efficacy_likelihood.R")
      
      TRUE
    })
    
    parallel::clusterExport(cl, c("params_list", "metadata_df", "base_inputs", "gen_bs",
                                  "n_particles", "n_threads", "burnints", "threshold", "tstep",
                                  "t_liverstage", "country_to_run", "VB", "divide", 
                                  # 'grid_search_list', 
                                  'observed_efficacy'
                                  ),
    envir = environment())
    
    # Grid search 
    grid_search_outputs <- parallel::clusterApply(cl,
                                                  params_list, 
                                                  function(params_row){
                                                    startsim <- Sys.time()
                                                    o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                                                                               metadata_df,
                                                                               base_inputs,
                                                                               output_dir = 'R:/Kelly/src/fit_smc/simulation_outputs',
                                                                               return_parasitemia = TRUE,
                                                                               save_outputs = FALSE)
                                                    endsim <- Sys.time()
                                                    message('time to do sim: ', endsim - startsim)
                                                    message('finished simulation')
                                                    filt <- o$infection_records %>% filter(time_ext >=0) %>% # filtering out the bites that happened before the start of FU
                                                      # remove any infections that occurred within 7 days 
                                                      group_by(rid) %>%
                                                      arrange(rid, detection_day) %>%
                                                      mutate(previous_detday = lag(detection_day),
                                                             diff = detection_day - previous_detday) %>%
                                                      filter(diff > 7 | is.na(diff)) %>% select(-diff, -previous_detday)
                                                    
                                                    eff <- calc_smc_efficacy_cumul(filt,
                                                                                   params_row,
                                                                                   by_week = TRUE)
                                                    
                                                    eff$sim_id <- params_row$sim_id
                                                    eff$repnum <- params_row$repnum
                                                    
                                                    matched <- observed_efficacy %>%
                                                      mutate(weeks_since_smc = ceiling(day_since_smc / 7)) %>%
                                                      group_by(weeks_since_smc) %>%
                                                      summarize(observed_efficacy = mean(efficacy)) %>%
                                                      mutate(observed_efficacy = ifelse(observed_efficacy == 1, 0.999, observed_efficacy)) %>%
                                                      left_join(eff %>% select(weeks_since_smc, efficacy) %>%
                                                                  rename(predicted_efficacy = efficacy), 
                                                                by = 'weeks_since_smc') %>% ungroup()
                                                    
                                                    # Remove NAs before calculating likelihood
                                                    matched_complete <- matched %>%
                                                      filter(!is.na(observed_efficacy), !is.na(predicted_efficacy))
                                                    
                                                    # Calculate mean least squares
                                                    mls <- mean((matched_complete$observed_efficacy - matched_complete$predicted_efficacy)^2)
                                                    
                                                    return(list(mls = mls, 
                                                                efficacy = eff,
                                                                params = params_row))
                                                  }
    )
    # results2 <- parallel::clusterApply(cl,
    #                                    params_list,
    #                                    function(params_row) {
    #                                      o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
    #                                                                 metadata_df,
    #                                                                 base_inputs,
    #                                                                 output_dir = 'R:/Kelly/synergy_orderly/src/fit_smc/simulation_outputs/',
    #                                                                 # allow_superinfections = TRUE,
    #                                                                 return_parasitemia = TRUE,
    #                                                                 save_outputs = FALSE)
    #                                      message('finished simulation')
    #                                      eff <- calc_smc_efficacy(o$infection_records,
    #                                                               params_row,
    #                                                               by_week = TRUE)
    #                                      eff_daily <- calc_smc_efficacy(o$infection_records,
    #                                                                     params_row,
    #                                                                     by_week = FALSE)
    #                                      eff$sim_id <- params_row$sim_id
    #                                      eff_daily$sim_id <- params_row$sim_id
    #                                      
    #                                      # Efficacy with cumulative proportion
    #                                      eff_cumul <- calc_smc_efficacy_cumul(o$infection_records,
    #                                                               params_row,
    #                                                               by_week = TRUE)
    #                                      eff_daily_cumul <- calc_smc_efficacy_cumul(o$infection_records,
    #                                                                     params_row,
    #                                                                     by_week = FALSE)
    #                                      eff_cumul$sim_id <- params_row$sim_id
    #                                      eff_daily_cumul$sim_id <- params_row$sim_id
    # 
    #                                      return(list(efficacy_weekly = eff,
    #                                                  efficacy_daily = eff_daily,
    #                                                  efficacy_weekly_cumul = eff_cumul, 
    #                                                  efficacy_daily_cumul = eff_daily_cumul,
    #                                                  parasitemia = o$parasitemia_data %>% mutate(sim_id = params_row$sim_id),
    #                                                  infection_records = o$infection_records %>% mutate(sim_id = params_row$sim_id),
    #                                                  params = params_row))
    # 
    #                                    }
    # )
    parallel::stopCluster(cl)
  }
  # Save all results 
  saveRDS(grid_search_outputs, paste0("R:/Kelly/synergy_orderly/src/fit_smc/outputs/runs_best_smcpars",Sys.Date(),".rds"))
  # saveRDS(results2, paste0("R:/Kelly/synergy_orderly/src/fit_smc/outputs/test_fitted_params_smc_",Sys.Date(),".rds"))
  saveRDS(metadata_df, paste0('R:/Kelly/synergy_orderly/src/fit_smc/outputs/metadata_', Sys.Date(), '.rds'))
}
