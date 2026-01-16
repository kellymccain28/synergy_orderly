run_grid_rtss <- function(path = "R:/Kelly/synergy_orderly",
                          n_param_sets,
                          threshold = 5000,
                          treatment_prob = 0.9, # default is 1 (which gives children prophylaxis)
                          N = 1200){
  
  # Script to fit SMC parameters to Hayley's curve 
  library(lhs)
  library(odin2)
  library(dust2)
  library(purrr)
  library(dplyr)
  library(mgcv)
  library(umbrella)
  # library(orderly2)
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
  source(paste0(path, "/src/fit_smc/calculate_efficacy_likelihood.R"))
  
  
  trial_ts = 365+90#*3# trial timesteps in cohort simulation (inte)
  # sim_allow_superinfections = TRUE # TRUE or FALSE
  country_to_run = 'generic'
  country_short = 'g'
  n_param_sets = n_param_sets
  N = N
  vax_day = -1 # unlike the model sim, this is in days (not timesteps)
  
  n_particles = 1L
  n_threads = 1L
  burnints = 50
  # threshold = 1000
  tstep = 1
  t_liverstage = 0 #7 # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC267587/
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
    country_short = country_short,
    treatment_probability = treatment_probability, 
    successful_treatment_probability = successful_treatment_probability
  )
  # # Generate LHS samples
  # Set up grid of parameter ranges
  # param_ranges <- list(
  #   alpha_ab = c(1.3, 1.9),
  #   beta_ab = c(2, 4),
  #   vmin = c(0, 0.01)
  # )
  # A <- randomLHS(n_param_sets, 3)
  # # Scale to parameter ranges
  # params_df <- data.frame(
  #   alpha_ab = qunif(A[,1], param_ranges$alpha_ab[1], param_ranges$alpha_ab[2]),
  #   beta_ab = qunif(A[,2], param_ranges$beta_ab[1], param_ranges$beta_ab[2]),
  #   vmin = qunif(A[,3], param_ranges$vmin[1], param_ranges$vmin[2])
  # )
  
  # params_df <- data.frame(
  #   alpha_ab = c(rep(1.476746, n_param_sets/3), rep(1.38, n_param_sets/3), rep(1.51, n_param_sets/3)),#1.285119 qunif(A[,1], param_ranges$alpha_ab[1], param_ranges$alpha_ab[2]),
  #   beta_ab = c(rep(5.83, n_param_sets/3), rep(5.83, n_param_sets/3), rep(6, n_param_sets/3)),#2.925123 qunif(A[,2], param_ranges$beta_ab[1], param_ranges$beta_ab[2]),
  #   vmin = 0#0.035372236 qunif(A[,3], param_ranges$vmin[1], param_ranges$vmin[2])
  # )
  # params_df <- data.frame(
  #   alpha_ab = rep(1.59, n_param_sets),#c(rep(1.83, n_param_sets/2),
  #   beta_ab = rep(2.46, n_param_sets),#c(rep(4.18, n_param_sets/2),
  #   vmin = rep(0.351, n_param_sets)#c(rep(0.322, n_param_sets/2)
  # )
  params_df <- rbind(
    # from post changing to mu parameterization 13/1
    # data.frame(alpha_ab = 1.73, beta_ab = 2.37, vmin = 0.000122),
    data.frame(alpha_ab = 1.77, beta_ab = 2.63, vmin = 0.000513)# this is the best one (see outputs from 01-13_4 and 01-14)
    # from pre changing to mu parameterization from prob (1-p) 13/1
    # data.frame(alpha_ab = 1.65, beta_ab = 3.26, vmin = 0.00395),
    # data.frame(alpha_ab = 1.48, beta_ab = 2.46, vmin = 0.00225)
  )
  params_df <- params_df %>%
    slice(rep(row_number(), each = n_param_sets))
  params_df <- params_df %>%
    mutate(
    max_SMC_kill_rate = rep(0, n_param_sets),
    lambda = rep(0, n_param_sets),
    kappa = rep(0, n_param_sets)
  )
  params_df$sim_id <- paste0('parameter_set_', rownames(params_df),"_", country_to_run, "_", treatment_probability)
  
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
  # saveRDS(parameters_df, 'parameters_df.rds')
  # Make metadata
  # Day of intervention (0 = start of follow-up; - values are before follow-up; + values after follow-up) - where 0 is also end of burnin
  # these get converted later to the correct directon - i.e. vaccine before follow-up will be +, smc before follow up will be -
  # vax_day is the 3rd primary dose (when we assume that efficacy begins)
  
  metadata_df <- data.frame(
    rid = 1:N,
    vaccination_day = vax_day,#sample(-1:0, N, replace = TRUE), the model takes abs(days before start of follow-up)
    PEV = c(rep(1, N/2), rep(0, N/2)),
    SMC = 0
  ) %>%
    mutate(arm = case_when(
      PEV == 1 & SMC == 1 ~ 'both',
      PEV == 1 & SMC == 0 ~ 'rtss',
      PEV == 0 & SMC == 1 ~ 'smc',
      TRUE ~ 'none')) %>%
    mutate(t_to_boost1 = 365,
           t_to_boost2 = 730,
           country = country_to_run) %>%
    mutate(rid_original = paste0(country_short, sprintf("%04d", rid)),
           country = 'generic',
           v1_date = as.Date('2017-04-01'))
  
  # "Observed" efficacy from White model 
  observed_efficacy_rtss <- readRDS(paste0(path, '/src/fit_rtss/observed_rtss_efficacy_months.rds'))
  
  # Run simulation
  cluster_cores <- Sys.getenv("CCP_NUMCPUS")
  # if (cluster_cores == "") {
  #   cluster_cores <- 8
  # }
  
  if (cluster_cores == "") {
    message("running in serial (on a laptop?)")
    message("Running ", nrow(params_df), " simulations sequentially")
    

    results2 <- lapply(params_list,
                       function(params_row) {
                         o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                                                    metadata_df,
                                                    base_inputs,
                                                    output_dir = 'R:/Kelly/src/fit_rtss/outputs',
                                                    return_parasitemia = FALSE,
                                                    save_outputs = FALSE)
                         message('finished simulation')
                         o$infection_records$sim_id <- params_row$sim_id
                         
                         infs <- o$infection_records %>%
                           # filter so that the follow-up time starts from 21 days post-vaccination 
                           filter(detection_day - vax_day >= 21) %>%
                           # remove any infections that occurred within 7 days 
                           group_by(rid) %>%
                           arrange(rid, detection_day) %>%
                           mutate(previous_detday = lag(detection_day),
                                  diff = detection_day - previous_detday) %>%
                           filter(diff > 7 | is.na(diff)) %>% select(-diff, -previous_detday)
                         
                         eff <- calc_rtss_efficacy(infs)
                         
                         matched <- observed_efficacy_rtss %>%
                           left_join(eff %>% select(months_since_rtss, efficacy) %>%
                                       rename(predicted_efficacy = efficacy), 
                                     by = 'months_since_rtss')
                         
                         # Remove NAs before calculating likelihood
                         matched_complete <- matched %>%
                           filter(!is.na(observed_efficacy), !is.na(predicted_efficacy))
                         
                         # Calculate mean least squares
                         mls <- mean((matched_complete$observed_efficacy - matched_complete$predicted_efficacy)^2, na.rm = TRUE)
                         
                         return(list(infection_records = o$infection_records,
                                     efficacy_weekly = eff,
                                     mls = mls, 
                                     params = params_row %>% select(-p_bite)))
                       }
    )
    
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
      # library(orderly2)
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
      source("R:/Kelly/synergy_orderly/src/fit_smc/calculate_efficacy_likelihood.R")
      
      TRUE
    })
    
    parallel::clusterExport(cl, c("params_list", "metadata_df", "base_inputs", "gen_bs",
                                  "n_particles", "n_threads", "burnints", "threshold", "tstep",
                                  "t_liverstage", "country_to_run", "VB", "divide",
                                  "observed_efficacy_rtss"
                                  ),
                            envir = environment())
    
    # message('starting optimization')

    results2 <- parallel::clusterApply(cl,
                                       params_list,
                                       function(params_row) {
                                         o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                                                                    metadata_df,
                                                                    base_inputs,
                                                                    output_dir = 'R:/Kelly/src/fit_rtss/outputs',
                                                                    return_parasitemia = FALSE,
                                                                    save_outputs = FALSE)
                                         message('finished simulation')
                                         o$infection_records$sim_id <- params_row$sim_id
                                         
                                         infs <- o$infection_records %>%
                                           # filter so that the follow-up time starts from 21 days post-vaccination 
                                           filter(detection_day - vax_day >= 21) %>%
                                           # remove any infections that occurred within 7 days 
                                           group_by(rid) %>%
                                           arrange(rid, detection_day) %>%
                                           mutate(previous_detday = lag(detection_day),
                                                  diff = detection_day - previous_detday) %>%
                                           filter(diff > 7 | is.na(diff)) %>% select(-diff, -previous_detday)
                                         
                                         eff <- calc_rtss_efficacy(infs)
                                         
                                         matched <- observed_efficacy_rtss %>%
                                           left_join(eff %>% select(months_since_rtss, efficacy) %>%
                                                       rename(predicted_efficacy = efficacy), 
                                                     by = 'months_since_rtss')
                                         
                                        # Remove NAs before calculating likelihood
                                         matched_complete <- matched %>%
                                           filter(!is.na(observed_efficacy), !is.na(predicted_efficacy))
                                         
                                         # Calculate mean least squares
                                         mls <- mean((matched_complete$observed_efficacy - matched_complete$predicted_efficacy)^2, na.rm = TRUE)
                                         
                                         return(list(infection_records = o$infection_records,
                                                     efficacy_weekly = eff,
                                                     mls = mls, 
                                                     params = params_row %>% select(-p_bite)))
                                       }
    )
    parallel::stopCluster(cl)
  }
  
  output_dir = paste0(path, '/src/fit_rtss/outputs/outputs_', Sys.Date())
  # If base directory doesn't exist, create it
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  } else {
    # If it exists, find an available numbered version
    counter <- 2
    new_dir <- paste0(output_dir, "_", counter)
    while (dir.exists(new_dir)) {
      counter <- counter + 1
      new_dir <- paste0(output_dir, "_", counter)
    }
    output_dir <- new_dir
    dir.create(output_dir, recursive = TRUE)
  }
  # saveRDS(results2, paste0(path, '/src/fit_rtss/outputs/rtss_grid_',Sys.Date(),'.rds'))
  # infectionrecords <- purrr::map_df(results2, "infection_records")
  efficacy <- purrr::map_df(results2, 'efficacy_weekly')
  params <- purrr::map_df(results2, 'params')
  mls <- lapply(results2, function(x) x$mls)
  mlsunlist <- unlist(mls)
  best <- which(mlsunlist == min(mlsunlist, na.rm = TRUE))[1]
  infectionrecords <- results2[[best]]$infection_records
  
  saveRDS(base_inputs, paste0(output_dir, '/base_inputs.rds'))
  saveRDS(params, paste0(output_dir, '/parameters.rds'))
  saveRDS(infectionrecords, paste0(output_dir, '/infectionrecords_rtss.rds'))
  saveRDS(efficacy, paste0(output_dir, '/efficacy_rtss.rds'))
  
  message(output_dir)
  # saveRDS(mls, paste0(output_dir, '/mls.rds'))
}


