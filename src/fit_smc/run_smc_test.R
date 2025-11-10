# script to quickly test the smc parameters from optimization after running fit_smc.R

run_smc_test <- function(path = "R:/Kelly/synergy_orderly",
                        n_param_sets,
                        N = 1200){
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
  source(paste0(path, "/src/sim_cohort_grid/cohort_sim_utils.R"))
  # Source processing functions
  source(paste0(path, "/shared/likelihood.R"))
  
  trial_ts = 80# trial timesteps in cohort simulation (inte)
  sim_allow_superinfections = TRUE # TRUE or FALSE
  country_to_run = 'generic'# BF or Mali, or if generic, then 'generic' which means that metadata_df is different.
  country_short = 'g'
  n_param_sets = n_param_sets
  N = N
  vax_day = 0 # unlike hte model sim, this is in days (not timesteps)
  
  n_particles = 1L
  n_threads = 1L
  burnints = 30
  threshold = 5000
  tstep = 1
  t_liverstage = 8
  VB = 1e6
  divide = if(tstep == 1) 2 else 1
  
  # Set up base inputs (these don't vary across parameter sweep)
  base_inputs <- list(
    trial_timesteps = trial_ts,
    burnin = burnints,
    threshold = threshold,
    VB = VB,
    tstep = tstep,
    t_liverstage = t_liverstage,
    country = country_to_run
  )
  
  # Set up grid of parameters
  param_ranges <- list(
    max_SMC_kill_rate = c(1, 10),# parasites per uL per 2-day timestep
    lambda = c(5, 50),
    kappa = c(0.1, 9)
  )
  # Generate LHS samples
  A <- randomLHS(n_param_sets, 3)
  # Scale to parameter ranges
  params_df <- data.frame(
    max_SMC_kill_rate = qunif(A[,1], param_ranges$max_SMC_kill_rate[1], param_ranges$max_SMC_kill_rate[2]),
    lambda = qunif(A[,2], param_ranges$lambda[1], param_ranges$lambda[2]),
    kappa = qunif(A[,3], param_ranges$kappa[1], param_ranges$kappa[2])#,
    # season_start_day = round(qunif(A[,4], param_ranges$season_start_day[1], param_ranges$season_start_day[2]))
  )
  # "fitted" parameter values for SMC
  # params_df <- params_df <- data.frame(
  #   max_SMC_kill_rate = rep(9, n_param_sets),
  #   lambda = rep(13.04, n_param_sets),
  #   kappa = rep(9, n_param_sets)
  # )
  params_df$sim_id <- paste0('parameter_set_', rownames(params_df),"_", country_to_run, "_", sim_allow_superinfections)
  
  prob_bite_generic <- readRDS(paste0(path, '/archive/fit_rainfall/20251009-144330-1d355186/prob_bite_generic.rds'))
  prob_bite_generic$prob_infectious_bite = 0.3
  p_bitevector <- calc_lagged_vectors(prob_bite_generic, 0, burnints = burnints) # no lagged values
  
  params_df$lag_p_bite <- 0
  params_df$season_start_day <- 0
  # SMC delivery
  params_df <- params_df %>%
    rowwise() %>%
    mutate(smc_dose_days = 10#list(c(seq(season_start_day, season_start_day + 120 - 1, 50),
           # seq(season_start_day + 365, season_start_day + 365 + 120 - 1, 50),
           # seq(season_start_day + 365*2, season_start_day  + 365*2 + 120 - 1,50)))
    ) %>%
    ungroup()
  
  parameters_df <- params_df %>%
    mutate(
      p_bite = purrr::map(lag_p_bite, ~p_bitevector[[paste0("lag_", .x)]])
    )
  
  # Make list of parameters instead of df
  params_list <- split(parameters_df, seq(nrow(parameters_df))) #%>%
  saveRDS(parameters_df, 'parameters_df.rds')
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
  
  
  # Run simulation
  cluster_cores <- Sys.getenv("CCP_NUMCPUS")
  if (cluster_cores == "") {
    cluster_cores <- 8
  }
  
  if (cluster_cores == "") {
    message("running in serial (on a laptop?)")
    message("Running ", nrow(params_df), " simulations sequentially")
    
    results2 <- lapply(params_list,
                       function(params_row){
                         o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                                                    metadata_df,
                                                    base_inputs,
                                                    output_dir = 'R:/Kelly/src/fit_smc/simulation_outputs',
                                                    allow_superinfections = TRUE,
                                                    return_parasitemia = FALSE,
                                                    save_outputs = FALSE)

                         eff <- calc_smc_efficacy(o$infection_records,
                                                  params_row,
                                                  by_week = TRUE)
                         eff_daily <- calc_smc_efficacy(o$infection_records,
                                                        params_row,
                                                        by_week = FALSE)
                         eff$sim_id <- params_row$sim_id
                         eff_daily$sim_id <- params_row$sim_id

                         return(list(efficacy_weekly = eff,
                                     efficacy_daily = eff_daily,
                                     params = params_row))
                       })
    
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
      
      source('R:/Kelly/synergy_orderly/src/sim_cohort_grid/cohort_sim_utils.R')
      source('R:/Kelly/synergy_orderly/shared/helper_functions.R')
      source("R:/Kelly/synergy_orderly/shared/rtss.R")
      source("R:/Kelly/synergy_orderly/shared/likelihood.R")
      source("R:/Kelly/synergy_orderly/src/fit_smc/calculate_efficacy_likelihood.R")
      
      TRUE
    })
    
    parallel::clusterExport(cl, c("params_list", "metadata_df", "base_inputs", "gen_bs",
                                  "n_particles", "n_threads", "burnints", "threshold", "tstep",
                                  "t_liverstage", "country_to_run", "VB", "divide"
    ),
    envir = environment())
    
    results2 <- parallel::clusterApply(cl,
                                       params_list,
                                       function(params_row) {
                                         o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                                                                    metadata_df,
                                                                    base_inputs,
                                                                    output_dir = 'R:/Kelly/src/fit_smc/simulation_outputs',
                                                                    allow_superinfections = TRUE,
                                                                    return_parasitemia = TRUE,
                                                                    save_outputs = FALSE)
                                         message('finished simulation')
                                         eff <- calc_smc_efficacy(o$infection_records,
                                                                  params_row,
                                                                  by_week = TRUE)
                                         eff_daily <- calc_smc_efficacy(o$infection_records,
                                                                        params_row,
                                                                        by_week = FALSE)
                                         eff$sim_id <- params_row$sim_id
                                         eff_daily$sim_id <- params_row$sim_id

                                         return(list(efficacy_weekly = eff,
                                                     efficacy_daily = eff_daily,
                                                     parasitemia = o$parasitemia_data %>% mutate(sim_id = params_row$sim_id),
                                                     infection_records = o$infection_records %>% mutate(sim_id = params_row$sim_id),
                                                     params = params_row))

                                       }
    )
    parallel::stopCluster(cl)
    
    # Save all results 
    saveRDS(results2, paste0("R:/Kelly/synergy_orderly/src/fit_smc/outputs/test_fitted_params_smc",Sys.Date(),".rds"))
  }
}

