run_fit_rtss <- function(path = "R:/Kelly/synergy_orderly",
                        n_param_sets,
                        treatment_prob = 0.9, # default is 1 (which gives children prophylaxis)
                        N = 1200){
  
  # Script to fit smc parameters to Hayley's curve 
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
  
  
  trial_ts = 365+100# trial timesteps in cohort simulation (inte)
  # sim_allow_superinfections = TRUE # TRUE or FALSE
  country_to_run = 'generic'
  country_short = 'g'
  n_param_sets = n_param_sets
  N = N
  vax_day = -1 # unlike the model sim, this is in days (not timesteps)
  
  n_particles = 1L
  n_threads = 1L
  burnints = 50
  threshold = 1000
  tstep = 1
  t_liverstage = 0#7 # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC267587/8
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
  
  params_df <- data.frame(
    max_SMC_kill_rate = rep(0, n_param_sets),
    lambda = rep(0, n_param_sets),
    kappa = rep(0, n_param_sets),
    alpha_ab = 1.38,#1.32, 
    beta_ab = 5.83#6.62
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
  saveRDS(parameters_df, 'parameters_df.rds')
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
  
  best_lhs <- data.frame(
    sim_id = 'parameter_set_1',
    alpha_ab = 1.38,
    beta_ab = 5.83,
    lag_p_bite = 0
  )
  # best_lhs <- pars[pars$sim_id %in% top_runs$sim_id,]
  best_lhs <- best_lhs %>%
    mutate(
      p_bite = purrr::map(lag_p_bite, ~p_bitevector[[paste0("lag_", .x)]]),
      smc_dose_days = 10
    )
  best_lhs_list <- split(best_lhs, seq(nrow(best_lhs)))
  
  # Run simulation
  cluster_cores <- Sys.getenv("CCP_NUMCPUS")
  if (cluster_cores == "") {
    cluster_cores <- 8
  }
  
  if (cluster_cores == "") {
    message("running in serial (on a laptop?)")
    message("Running ", nrow(params_df), " simulations sequentially")
    
    # For optimization
    optim_results <- lapply(
      best_lhs_list,
      function(start){
        initial_params <- c(start$alpha_ab,
                            start$beta_ab)
        lower_bounds <- c(0.5, 4) # alpha, beta
        upper_bounds <- c(3, 8)
        
        # Track evaluations
        n_evals <- 0
        eval_history <- list()
        
        objective <- function(params) {
          tryCatch({
            n_evals <<- n_evals + 1
            
            message(sprintf("\n=== Evaluation %d ===", n_evals))
            message(sprintf("Params: alpha=%.4f, beta=%.4f", 
                            params[1], params[2]))
            
            params_tibble <- data.frame(
              max_SMC_kill_rate = 0,
              lambda = 0,
              kappa = 0,
              alpha_ab = params[1],
              beta_ab = params[2],
              lag_p_bite = 0,
              smc_dose_days = start$smc_dose_days,
              sim_id = start$sim_id
            )
            params_tibble$p_bite <- list(start$p_bite)
            
            mls <- calculate_efficacy_likelihood_rtss(params_tibble,
                                                      metadata_df,
                                                      base_inputs,
                                                      observed_efficacy_rtss )
            
            message('Evaluation ', n_evals, ': mls = ', round(mls, 4))
            
            # Store history with tibble
            eval_history[[n_evals]] <<- list(
              params_tibble = params_tibble[4:5],
              mls = mls
            )
            
            return(mls)
          }, error = function(e) {
            message("Error in objective function: ", e$message)
            print(traceback())
            stop(e)
          })
        }
        
        # Run optimization with STRICT iteration limit
        fit <- optim(
          par = initial_params,
          fn = objective,
          method = "L-BFGS-B",
          lower = lower_bounds,
          upper = upper_bounds,
          control = list(
            maxit = 100,  # Hard limit
            trace = 1,
            factr = 1e6  # Loose convergence 
          )
        )
        
        return(list(
          starting_point_id = start$sim_id,
          initial_params = initial_params,
          final_params = fit$par,
          mls = fit$value,
          convergence = fit$convergence,
          n_evaluations = n_evals,
          eval_history = eval_history,
          fit = fit
        ))
        
      })
    # results2 <- lapply(params_list,
    #                    function(params_row){
    #                      o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
    #                                                 metadata_df,
    #                                                 base_inputs,
    #                                                 output_dir = 'R:/Kelly/src/fit_rtss/outputs',
    #                                                 # allow_superinfections = TRUE,
    #                                                 return_parasitemia = FALSE,
    #                                                 save_outputs = FALSE)
    #                      message('finished simulation')
    #                      o$infection_records$sim_id <- params_row$sim_id
    #                      
    #                      eff <- calc_rtss_efficacy(o$infection_records)
    #                      eff_cumul <- calc_rtss_efficacy_cumul(o$infection_records)
    #                      
    #                      return(list(infection_records = o$infection_records, 
    #                                  efficacy_weekly = eff,
    #                                  efficacy_weekly_cumul = eff_cumul,
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
                                  "observed_efficacy_rtss", "best_lhs_list"),
                            envir = environment())
    
    message('starting optimization')
    # For optimization
    # For optimization
    optim_results <- parallel::clusterApply(cl,
                                            best_lhs_list,
                                            function(start){
                                              initial_params <- c(start$alpha_ab,
                                                                  start$beta_ab)
                                              lower_bounds <- c(1, 3) # alpha, beta
                                              upper_bounds <- c(4, 7)
                                              
                                              # Track evaluations
                                              n_evals <- 0
                                              eval_history <- list()
                                              
                                              objective <- function(params) {
                                                tryCatch({
                                                  n_evals <<- n_evals + 1
                                                  
                                                  message(sprintf("\n=== Evaluation %d ===", n_evals))
                                                  message(sprintf("Params: alpha=%.4f, beta=%.4f", 
                                                                  params[1], params[2]))
                                                  
                                                  params_tibble <- data.frame(
                                                    max_SMC_kill_rate = 0,
                                                    lambda = 0,
                                                    kappa = 0,
                                                    alpha_ab = params[1],
                                                    beta_ab = params[2],
                                                    lag_p_bite = 0,
                                                    smc_dose_days = start$smc_dose_days,
                                                    sim_id = start$sim_id
                                                  )
                                                  params_tibble$p_bite <- list(start$p_bite)
                                                  
                                                  mls <- calculate_efficacy_likelihood_rtss(params_tibble,
                                                                                            metadata_df,
                                                                                            base_inputs,
                                                                                            observed_efficacy_rtss )
                                                  
                                                  message('Evaluation ', n_evals, ': mls = ', round(mls, 4))
                                                  
                                                  # Store history with tibble
                                                  eval_history[[n_evals]] <<- list(
                                                    params_tibble = params_tibble[4:5],
                                                    mls = mls
                                                  )
                                                  
                                                  return(mls)
                                                }, error = function(e) {
                                                  message("Error in objective function: ", e$message)
                                                  print(traceback())
                                                  stop(e)
                                                })
                                              }
                                              
                                              # Run optimization with STRICT iteration limit
                                              fit <- optim(
                                                par = initial_params,
                                                fn = objective,
                                                method = "L-BFGS-B",
                                                lower = lower_bounds,
                                                upper = upper_bounds,
                                                control = list(
                                                  maxit = 100,  # Hard limit
                                                  trace = 1,
                                                  factr = 1e7  # Loose convergence 
                                                )
                                              )
                                              
                                              return(list(
                                                starting_point_id = start$sim_id,
                                                initial_params = initial_params,
                                                final_params = fit$par,
                                                mls = fit$value,
                                                convergence = fit$convergence,
                                                n_evaluations = n_evals,
                                                eval_history = eval_history,
                                                fit = fit
                                              ))
                                              
                                            })
    # results2 <- parallel::clusterApply(cl,
    #                                    params_list,
    #                                    function(params_row) {
    #                                      o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
    #                                                                 metadata_df,
    #                                                                 base_inputs,
    #                                                                 output_dir = 'R:/Kelly/src/fit_rtss/outputs',
    #                                                                 # allow_superinfections = TRUE,
    #                                                                 return_parasitemia = FALSE,
    #                                                                 save_outputs = FALSE)
    #                                      message('finished simulation')
    #                                      o$infection_records$sim_id <- params_row$sim_id
    #                                      
    #                                      eff <- calc_rtss_efficacy(o$infection_records)
    #                                      eff_cumul <- calc_rtss_efficacy_cumul(o$infection_records)
    #                                      
    #                                      return(list(infection_records = o$infection_records, 
    #                                                  efficacy_weekly = eff,
    #                                                  efficacy_weekly_cumul = eff_cumul,
    #                                                  params = params_row))
    #                                    }
    # )
    parallel::stopCluster(cl)
  }
  
  saveRDS(optim_results, paste0(path, '/src/fit_rtss/outputs/optimization_results_2111.rds'))
  # infectionrecords <- purrr::map_df(results2, "infection_records")
  # efficacy <- purrr::map_df(results2, 'efficacy_weekly')
  # efficacy_cumul <- purrr::map_df(results2, 'efficacy_weekly_cumul')
  # params <- purrr::map_df(results2, 'parameters')
  # 
  # saveRDS(params, paste0(path, '/src/fit_rtss/outputs/parameters_', Sys.Date(), '.rds'))
  # saveRDS(infectionrecords, paste0(path, '/src/fit_rtss/outputs/infectionrecords_rtss_', Sys.Date(), '.rds'))
  # saveRDS(efficacy, paste0(path, '/src/fit_rtss/outputs/efficacy_rtss_', Sys.Date(), '.rds'))
  # saveRDS(efficacy_cumul, paste0(path, '/src/fit_rtss/outputs/efficacy_rtss_cumul_', Sys.Date(), '.rds'))
}


