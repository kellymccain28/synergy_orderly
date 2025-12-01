run_fit_smc <- function(path = "R:/Kelly/synergy_orderly",
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
  source(paste0(path, "/shared/likelihood.R"))
  source(paste0(path, "/src/fit_smc/calculate_efficacy_likelihood.R"))
  
  trial_ts = 80# trial timesteps in cohort simulation (inte)
  # sim_allow_superinfections = TRUE # TRUE or FALSE
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
  
  # Set up grid of parameters
  param_ranges <- list(
    max_SMC_kill_rate = c(1, 10),# parasites per uL per 2-day timestep
    lambda = c(5, 50),
    kappa = c(0.01, 5)
  )
  # # Generate LHS samples
  # A <- randomLHS(n_param_sets, 3)
  # # Scale to parameter ranges
  # params_df <- data.frame(
  #   max_SMC_kill_rate = qunif(A[,1], param_ranges$max_SMC_kill_rate[1], param_ranges$max_SMC_kill_rate[2]),
  #   lambda = qunif(A[,2], param_ranges$lambda[1], param_ranges$lambda[2]),
  #   kappa = qunif(A[,3], param_ranges$kappa[1], param_ranges$kappa[2])
  # )
  # "fitted" parameter values for SMC
  params_df <- params_df <- data.frame(
    max_SMC_kill_rate = rep(2.89, n_param_sets),#rep(3, n_param_sets),
    lambda =rep(17, n_param_sets),# rep(13.08, n_param_sets),
    kappa = rep(0.28, n_param_sets)#rep(0.43, n_param_sets)
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
  
  
  # Get best LHS values from previous runs of LHS -- 
  observed_efficacy <- read.csv(paste0(path, '/shared/smc_fits_hayley.csv')) %>%
    mutate(weeks_since_smc = ceiling(day_since_smc / 7)) %>%
    group_by(weeks_since_smc) %>%
    mutate(efficacy = mean(efficacy)) # because teh first week should be 1 not 0
  # `efficacy_weekly_2025-10-15` <- readRDS("R:/Kelly/synergy_orderly/src/fit_smc/outputs/efficacy_weekly_2025-10-15.rds") %>%
  #   mutate(sim_id = paste0(sim_id, '1015'))
  # `efficacy_weekly_2025-10-14` <- readRDS("R:/Kelly/synergy_orderly/src/fit_smc/outputs/efficacy_weekly_2025-10-14.rds") %>%
  #   mutate(sim_id = paste0(sim_id, '1014'))
  # `efficacy_weekly_2025-10-14_100` <- readRDS("R:/Kelly/synergy_orderly/src/fit_smc/outputs/efficacy_weekly_2025-10-14_100.rds") %>%
  #   mutate(sim_id = paste0(sim_id, '1014_1'))
  # efficacy_results <- bind_rows(`efficacy_weekly_2025-10-14`, `efficacy_weekly_2025-10-14_100`, `efficacy_weekly_2025-10-15`)
  # `parameters_2025-10-15` <- readRDS("R:/Kelly/synergy_orderly/src/fit_smc/outputs/parameters_2025-10-15.rds")%>%
  #   mutate(sim_id = paste0(sim_id, '1015'))
  # `parameters_2025-10-14` <- readRDS("R:/Kelly/synergy_orderly/src/fit_smc/outputs/parameters_2025-10-14.rds") %>%
  #   mutate(sim_id = paste0(sim_id, '1014'))
  # `parameters_2025-10-14_100` <- readRDS("R:/Kelly/synergy_orderly/src/fit_smc/outputs/parameters_2025-10-14_100.rds") %>%
  #   mutate(sim_id = paste0(sim_id, '1014_1'))
  # pars <- bind_rows(`parameters_2025-10-14`, `parameters_2025-10-14_100`, `parameters_2025-10-15`)

  # # Assume observed data has Normal measurement error
  # calculate_likelihood <- function(observed, predicted, sigma = 0.1) {
  #   # neg log-likelihood
  #   -sum(dnorm(observed, mean = predicted, sd = sigma, log = TRUE), na.rm = TRUE)
  # }
  # 
  # likelihoods <- efficacy_results %>%
  #   ungroup() %>%
  #   rename(predicted_efficacy = efficacy) %>%
  #   left_join(observed_efficacy, by = "weeks_since_smc") %>%
  #   group_by(sim_id) %>%
  #   summarise(
  #     neg_log_lik = calculate_likelihood(observed_efficacy,
  #                                        predicted_efficacy,
  #                                        sigma = 0.1),  # Adjust sigma based on your uncertainty
  #     .groups = 'drop'
  #   ) %>%
  #   mutate(
  #     likelihood = exp(-neg_log_lik),
  #     weight = likelihood / sum(likelihood)  # Normalize to get weights
  #   ) %>%
  #   arrange(neg_log_lik)
  # 
  # # Get top 10 runs
  # top_runs <- likelihoods %>% slice_head(n = 5)

  # efficacy_results %>%
  #   filter(sim_id %in% top_runs$sim_id &
  #            weeks_since_smc < 10
  #   ) %>%
  #   ggplot(aes(x = weeks_since_smc*7, y = efficacy, group = sim_id, color = as.factor(sim_id))) +
  #   geom_line(alpha = 0.5) +
  #   geom_point(data = observed_efficacy, aes(x = weeks_since_smc*7, y = observed_efficacy),
  #              color = "black", inherit.aes = FALSE) +
  #   geom_line(data = observed_efficacy, aes(x = weeks_since_smc*7, y = observed_efficacy),
  #             color = "black", inherit.aes = FALSE) +
  #   labs(title = "Top 10 runs vs observed SMC efficacy",
  #        y = "SMC efficacy", x = "Days since SMC") +
  #   theme_minimal()  + theme(legend.position = 'none')
  best_lhs <- data.frame(
    max_SMC_kill_rate = c(2.89, 3),#c(3,3,4,4,3,3,4,4) ,
    lambda = c(17.3,15),# c(13.08,14,13.08,14,13.08,14,13.08,14),
    kappa = c(0.278,0.454),#c(0.454,0.454,0.454,0.454, 0.5,0.5, 0.5,0.5),
    sim_id = c(1,2),#seq(1,4),
    lag_p_bite = rep(0, 2)#0
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
  # if (cluster_cores == "") {
  #   cluster_cores <- 8
  # }
  
  if (cluster_cores == "") {
    message("running in serial (on a laptop?)")
    message("Running ", nrow(params_df), " simulations sequentially")
    
    # For optimization
    optim_results <- lapply(
      best_lhs_list,
      function(start){
        initial_params <- c(start$max_SMC_kill_rate,
                            start$lambda,
                            start$kappa)
        lower_bounds <- c(1, 10, 0.01) # max, lambda, kappa
        upper_bounds <- c(16, 25, 2)

        # Track evaluations
        n_evals <- 0
        eval_history <- list()

        objective <- function(params) {
          n_evals <<- n_evals + 1
          
          message(sprintf("\n=== Evaluation %d ===", n_evals))
          message(sprintf("Params: max=%.4f, lambda=%.4f, kappa=%.4f", 
                          params[1], params[2], params[3]))

          params_tibble <- data.frame(
            max_SMC_kill_rate = params[1],
            lambda = params[2],
            kappa = params[3],
            lag_p_bite = 0,
            smc_dose_days = start$smc_dose_days,
            sim_id = start$sim_id
          )
          params_tibble$p_bite <- list(start$p_bite)

          mls <- calculate_efficacy_likelihood(params_tibble,
                                        metadata_df,
                                        base_inputs,
                                        observed_efficacy )
          
          message('Evaluation ', n_evals, ': mls = ', round(mls, 4))

          # Store history with tibble
          eval_history[[n_evals]] <<- list(
            params_tibble = params_tibble[1:3],
            mls = mls
          )
          
          return(mls)
        }

        # Run optimization with STRICT iteration limit
        fit <- optim(
          par = initial_params,
          fn = objective,
          method = "L-BFGS-B",
          lower = lower_bounds,
          upper = upper_bounds,
          control = list(
            maxit = 800,  # Hard limit
            trace = 1,
            factr = 1e7  # Loose convergence 
          )
        )

        return(list(
          starting_point_id = start$sim_id,
          initial_params = initial_params,
          final_params = fit$par,
          mean_least_squares = fit$value,
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
    #                                                 output_dir = 'R:/Kelly/src/fit_smc/simulation_outputs',
    #                                                 allow_superinfections = TRUE,
    #                                                 return_parasitemia = FALSE,
    #                                                 save_outputs = FALSE)
    # 
    #                      eff <- calc_smc_efficacy(o$infection_records,
    #                                               params_row,
    #                                               by_week = TRUE)
    #                      eff_daily <- calc_smc_efficacy(o$infection_records,
    #                                                     params_row,
    #                                                     by_week = FALSE)
    #                      eff$sim_id <- params_row$sim_id
    #                      eff_daily$sim_id <- params_row$sim_id
    # 
    #                      return(list(efficacy_weekly = eff,
    #                                  efficacy_daily = eff_daily,
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
      source("R:/Kelly/synergy_orderly/shared/likelihood.R")
      source("R:/Kelly/synergy_orderly/src/fit_smc/calculate_efficacy_likelihood.R")
      
      TRUE
    })
    
    parallel::clusterExport(cl, c("params_list", "metadata_df", "base_inputs", "gen_bs",
                                  "n_particles", "n_threads", "burnints", "threshold", "tstep",
                                  "t_liverstage", "country_to_run", "VB", "divide",
                                  "observed_efficacy",
                                  "best_lhs_list",  "param_ranges"
    ),
                            envir = environment())
    # for optimization
    optim_results <- parallel::clusterApply(cl,
                                            best_lhs_list,
                                            function(start){
                                              initial_params <- c(start$max_SMC_kill_rate,
                                                                  start$lambda,
                                                                  start$kappa)
                                              lower_bounds <- c(1, 10, 0.01) # max, lambda, kappa
                                              upper_bounds <- c(16, 25, 2)
                                              
                                              # Track evaluations
                                              n_evals <- 0
                                              eval_history <- list()
                                              
                                              objective <- function(params) {
                                                n_evals <<- n_evals + 1
                                                
                                                message(sprintf("\n=== Evaluation %d ===", n_evals))
                                                message(sprintf("Params: max=%.4f, lambda=%.4f, kappa=%.4f", 
                                                                params[1], params[2], params[3]))
                                                
                                                params_tibble <- data.frame(
                                                  max_SMC_kill_rate = params[1],
                                                  lambda = params[2],
                                                  kappa = params[3],
                                                  lag_p_bite = 0,
                                                  smc_dose_days = start$smc_dose_days,
                                                  sim_id = start$sim_id
                                                )
                                                params_tibble$p_bite <- list(start$p_bite)
                                                
                                                mls <- calculate_efficacy_likelihood(params_tibble,
                                                                                     metadata_df,
                                                                                     base_inputs,
                                                                                     observed_efficacy )
                                                
                                                message('Evaluation ', n_evals, ': mls = ', round(mls, 4))
                                                
                                                # Store history with tibble
                                                eval_history[[n_evals]] <<- list(
                                                  params_tibble = params_tibble[1:3],
                                                  mls = mls
                                                )
                                                
                                                return(mls)
                                              }
                                              
                                              # Run optimization with STRICT iteration limit
                                              fit <- optim(
                                                par = initial_params,
                                                fn = objective,
                                                method = "L-BFGS-B",#"L-BFGS-B",
                                                lower = lower_bounds,
                                                upper = upper_bounds,
                                                control = list(
                                                  maxit = 800,  # Hard limit
                                                  trace = 1,
                                                  factr = 1e5,      # Tighter tolerance
                                                  pgtol = 1e-8      # Gradient tolerance
                                                )
                                              )
                                              
                                              return(list(
                                                starting_point_id = start$sim_id,
                                                initial_params = initial_params,
                                                final_params = fit$par,
                                                mean_least_squares = fit$value,
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
    #                                                                 output_dir = 'R:/Kelly/src/fit_smc/simulation_outputs',
    #                                                                 allow_superinfections = TRUE,
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
    #                                      return(list(efficacy_weekly = eff,
    #                                                  efficacy_daily = eff_daily,
    #                                                  parasitemia = o$parasitemia_data,
    #                                                  params = params_row))
    # 
    #                                    }
    # )
    parallel::stopCluster(cl)
    
    # Save all results 
    saveRDS(optim_results, paste0('R:/Kelly/synergy_orderly/src/fit_smc/outputs/optimization_results_',Sys.Date(), '.rds'))
    # saveRDS(results2, "R:/Kelly/synergy_orderly/src/fit_smc/outputs/test_fitted_params_smc.rds")
  }
  
}


# Observed efficacy extracted from figure in SI of Hayley's paper 10.1016/S2214-109X(22)00416-8
# observed_efficacy <- read.csv(paste0(path, '/shared/smc_fits_hayley.csv'))
# lhs_parameters <- readRDS("R:/Kelly/synergy_orderly/src/fit_smc/outputs/test_fitted_params_smc_0611.rds")#lhs_parameters_0411.rds
# effweekly <- purrr::map_df(lhs_parameters, 'efficacy_weekly')
# effdaily <- purrr::map_df(lhs_parameters, 'efficacy_daily')
# lhspars <- purrr::map_df(lhs_parameters, 'params')
# # # # Look at efficacy output
# ggplot(effdaily %>% filter(days_since_smc < 70)) +
#   geom_point(aes(x = days_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.1) +
#   geom_line(aes(x = days_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.2) +
#   ylim(c(-0.5, 1)) +
#   geom_line(data = observed_efficacy, aes(x = day_since_smc, y = efficacy)) +
#   geom_point(data = observed_efficacy, aes(x = day_since_smc, y = efficacy)) +
#   theme_minimal() +
#   theme(legend.position = 'none')
# 
# ggplot(effweekly %>% filter(weeks_since_smc < 10)) +
#   geom_point(aes(x = weeks_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.1) +
#   geom_line(aes(x = weeks_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.2) +
#   ylim(c(-0.5, 1)) +
#   geom_line(data = observed_efficacy, aes(x = day_since_smc/7, y = efficacy)) +
#   geom_point(data = observed_efficacy, aes(x = day_since_smc/7, y = efficacy)) +
#   scale_x_continuous(breaks = seq(0,9,1)) +
#   theme_minimal() +
#   theme(legend.position = 'none')
# #
# # # Assume observed data has Normal measurement error
# calculate_likelihood <- function(observed, predicted, sigma = 0.1) {
#   # neg log-likelihood
#   -sum(dnorm(observed, mean = predicted, sd = sigma, log = TRUE), na.rm = TRUE)
# }
# 
# eff <- effweekly
# eff <- effdaily
# likelihoods <- eff %>%
#   ungroup() %>%
#   mutate(#day_since_smc = round(weeks_since_smc * 7,0),
#     day_since_smc = days_since_smc
#   ) %>%
#   rename(predicted_efficacy = efficacy) %>%
#   left_join(observed_efficacy %>% rename(observed_efficacy = efficacy), by = "day_since_smc") %>%
#   group_by(sim_id) %>%
#   summarise(
#     neg_log_lik = calculate_likelihood(observed_efficacy,
#                                        predicted_efficacy,
#                                        sigma = 0.1),  # Adjust sigma based on your uncertainty
#     .groups = 'drop'
#   ) %>%
#   mutate(
#     likelihood = exp(-neg_log_lik),
#     weight = likelihood / sum(likelihood)  # Normalize to get weights
#   ) %>%
#   arrange(neg_log_lik)
# 
# # Get top 10 runs
# top_runs <- likelihoods %>% slice_head(n = 10)
# 
# # Plot them
# eff %>%
#   filter(sim_id %in% top_runs$sim_id &
#            # days_since_smc < 10*7
#          weeks_since_smc < 10
#   ) %>%
#   ggplot(aes(x = weeks_since_smc*7, y = efficacy, group = sim_id, color = as.factor(sim_id))) +
#   # ggplot(aes(x = days_since_smc, y = efficacy, group = sim_id, color = as.factor(sim_id))) +
#   geom_line(alpha = 0.5) +
#   geom_point(data = observed_efficacy, aes(x = day_since_smc, y = efficacy),
#              color = "black", inherit.aes = FALSE) +
#   geom_line(data = observed_efficacy, aes(x = day_since_smc, y = efficacy),
#             color = "black", inherit.aes = FALSE) +
#   labs(title = "Top 10 runs vs observed SMC efficacy",
#        y = "SMC efficacy", x = "Days since SMC") +
# theme_minimal() + ylim(c(-0.5, 1)) + theme(legend.position = 'none')
# #
# # # parameter_set_293_generic_TRUE is closest
# # # max_SMC_kill_rate lambda kappa season_start_day sim_id                         lag_p_bite
# # # <dbl>             <dbl>  <dbl> <dbl>            <chr>                          <dbl>
# # # 6.41              44.2   0.425 150              parameter_set_293_generic_TRUE 0
# #
# # top 10
# # max_SMC_kill_rate lambda kappa season_start_day sim_id                         lag_p_bite
# # <dbl>  <dbl> <dbl>            <dbl> <chr>                               <dbl>
# #   1             12.7   41.2  8.20                98 parameter_set_114_generic_TRUE          0
# # 2              8.80  36.8  9.55                51 parameter_set_133_generic_TRUE          0
# # 3              7.99   9.99 0.235               41 parameter_set_248_generic_TRUE          0
# # 4              6.41  44.2  0.425              150 parameter_set_293_generic_TRUE          0
# # 5             24.5   30.0  3.71                99 parameter_set_420_generic_TRUE          0
# # 6             20.0   36.1  7.71                89 parameter_set_545_generic_TRUE          0
# # 7             21.0   32.1  4.15                37 parameter_set_802_generic_TRUE          0
# # 8             16.0   33.5  8.39                74 parameter_set_867_generic_TRUE          0
# # 9             14.4   33.5  9.87                94 parameter_set_934_generic_TRUE          0
# # 10             22.7   31.0  9.05               123 parameter_set_993_generic_TRUE          0
# 
# eff %>%
#   filter(sim_id %in% top_runs$sim_id &
#            days_since_smc < 10*5
#          # weeks_since_smc < 10
#   ) %>%
#   ggplot() +
#   # geom_line(aes(x = weeks_since_smc, y = inci_none, group = sim_id, color = as.factor(sim_id)),
#   #           alpha = 0.8, linetype = 2) +
#   # geom_line(aes(x = weeks_since_smc, y = inci_smc, group = sim_id, color = as.factor(sim_id)),
#   #           alpha = 0.8) +
#   geom_line(aes(x = days_since_smc, y = inci_none, group = sim_id, color = as.factor(sim_id)),
#             alpha = 0.8, linetype = 2) +
#   geom_line(aes(x = days_since_smc, y = inci_smc, group = sim_id, color = as.factor(sim_id)),
#             alpha = 0.8) +
#   labs(y = "incidence", x = "days since SMC") +
#   theme_minimal() #+ ylim(c(0,0.7))
# 
# weibull_survival <- function(t, max_SMC_kill_rate, lambda, kappa){
#   max_SMC_kill_rate * exp(-(t/lambda)^kappa)
# }
# 
# t_seq = seq(0, 60, length.out = 60)
# pars <- `parameters_2025-10-14` %>% filter(sim_id %in% top_runs$sim_id)
# 
# plot_data <- pars %>%
#   mutate(numid = readr::parse_number(sim_id)) %>%
#   mutate(label = paste0(numid, "\n", round(max_SMC_kill_rate,2), ', ', round(lambda,2),', ', round(kappa,2))) %>%
#   rowwise() %>%
#   do({
#     data.frame(
#       id = .$sim_id,
#       label = .$label,
#       time = t_seq,
#       survival = weibull_survival(t_seq, .$max_SMC_kill_rate, .$lambda, .$kappa)
#     )
#   }) %>%
#   ungroup()
# 
# ggplot(plot_data, aes(x = time, y = survival)) +
#   geom_line(color = "steelblue", linewidth = 1) +
#   facet_wrap(~ label) +
#   labs(
#     title = "Weibull Survival Functions",
#     x = "Time",
#     y = "Survival Probability S(t)"
#   ) +
#   theme_minimal() +
#   theme(
#     strip.text = element_text(size = 11, face = "bold"),
#     plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
#     legend.position = 'none'
#   )
# 
# 
# 
# 
# # Optimization
# optimization_results <- readRDS("R:/Kelly/synergy_orderly/src/fit_smc/outputs/optimization_results_20251113.rds") # 2210 does allow superinfections
# 
# best_refined <- optimization_results[[which.max(
#   sapply(optimization_results, function(x) x$log_likelihood)
# )]]
# 
# best_refined$log_likelihood
# best_refined$starting_point_id
# best_refined$initial_params #9.0928452 13.0471435  0.4141795
# best_refined$final_params #13 nov 3.0093347 13.0480119  0.4539757 #9.0739477 15.0000000  0.4586005 # when allowing lambda to go <15, 9.1078524 13.0757521  0.4277123
# best_refined$convergence#0
# best_refined$n_evaluations #231
# pars_to_test <- data.frame(
#   max_SMC_kill_rate = rep(best_refined$final_params[1], 10),
#   lambda = rep(best_refined$final_params[2],10),
#   kappa = rep(best_refined$final_params[3],10),
#   sim_id = seq(1, 10),
#   lag_p_bite = rep(0,10),
#   season_start_day = rep(50,10),
#   smc_dose_days = rep(10,10),
#   p_bite = I(rep(params_list[[1]]$p_bite, 10))
# )
# pars_to_test <- split(pars_to_test, seq(nrow(pars_to_test)))
# results2 <- lapply(pars_to_test,
#                    function(params_row){
#                      o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
#                                                 metadata_df,
#                                                 base_inputs,
#                                                 output_dir = 'R:/Kelly/src/fit_smc/simulation_outputs',
#                                                 allow_superinfections = TRUE,
#                                                 return_parasitemia = TRUE,
#                                                 save_outputs = FALSE)
# 
#                      eff <- calc_smc_efficacy(o$infection_records,
#                                               params_row,
#                                               by_week = TRUE)
#                      eff_daily <- calc_smc_efficacy(o$infection_records,
#                                                     params_row,
#                                                     by_week = FALSE)
#                      eff$sim_id <- params_row$sim_id
#                      eff_daily$sim_id <- params_row$sim_id
# 
#                      return(list(efficacy_weekly = eff,
#                                  efficacy_daily = eff_daily,
#                                  params = params_row,
#                                  parasitemia = o$parasitemia_data,
#                                  infection_records_smc = o$infection_records))
#                    })
# eff_weekly <- purrr::map_df(results2, 'efficacy_weekly')
# eff_daily <- purrr::map_df(results2, 'efficacy_daily')
# ggplot(eff_daily %>% filter(days_since_smc < 70)) +
#   geom_point(aes(x = days_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.8) +
#   geom_line(aes(x = days_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.2) +
#   ylim(c(-0.5, 1)) +
#   geom_line(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), alpha = 0.2) +
#   geom_point(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), alpha = 0.2) +
#   theme_minimal() +
#   theme(legend.position = 'none')
# ggplot(eff_weekly %>% filter(weeks_since_smc < 10)) +
#   ylim(c(-0.5, 1)) +
#   geom_line(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), alpha = 0.8) +
#   geom_point(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), alpha = 0.8) +
#   geom_point(aes(x = weeks_since_smc*7, y = efficacy, group = sim_id, color = as.factor(sim_id)), alpha = 0.8) +
#   geom_line(aes(x = weeks_since_smc*7, y = efficacy, group = sim_id, color = as.factor(sim_id)), alpha = 0.8) +
#   theme_minimal() +
#   theme(legend.position = 'none')
# 
# ggplot(eff_weekly %>% filter(weeks_since_smc < 10)) +
#     geom_line(aes(x = weeks_since_smc, y = inci_none, group = sim_id, color = as.factor(sim_id)),
#               alpha = 0.8, linetype = 2) +
#     geom_line(aes(x = weeks_since_smc, y = inci_smc, group = sim_id, color = as.factor(sim_id)),
#               alpha = 0.8) +
#     labs(y = "incidence", x = "weeks since SMC") +
#     theme_minimal() +#+ ylim(c(0,0.7))
#     theme(legend.position = 'none')
# 
# ggplot(eff_daily %>% filter(days_since_smc < 70)) +
#   geom_line(aes(x = days_since_smc, y = inci_none, group = sim_id, color = as.factor(sim_id)),
#             alpha = 0.8, linetype = 2) +
#   geom_line(aes(x = days_since_smc, y = inci_smc, group = sim_id, color = as.factor(sim_id)),
#             alpha = 0.8) +
#   labs(y = "incidence", x = "days since SMC") +
#   theme_minimal() +#+ ylim(c(0,0.7))
#   theme(legend.position = 'none')
# 
