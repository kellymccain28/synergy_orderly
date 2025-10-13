run_fit_smc <- function(path = "R:/Kelly/synergy_orderly",
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
  source(paste0(path, "/src/sim_cohort_grid/cohort_sim_utils.R"))
  # SOurce processing functions
  source(paste0(path, "/shared/likelihood.R"))
  # source('shared/get_cox_efficacy.R')
  # source("shared/format_model_output.R")
  # source("shared/get_incidence.R")
  
  trial_ts = 365# trial timesteps in cohort simulation (inte)
  sim_allow_superinfections = TRUE # TRUE or FALSE
  country_to_run = 'generic'# BF or Mali, or if generic, then 'generic' which means that metadata_df is different.
  country_short = 'g'
  n_param_sets = n_param_sets
  N = N
  vax_day = -10 # unlike hte model sim, this is in days (not timesteps)
  
  n_particles = 1L
  n_threads = 1L
  burnints = 50
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
    max_SMC_kill_rate = c(3, 25),# parasites per uL per 2-day timestep
    lambda = c(10, 50),
    kappa = c(0.1, 10)#,
    # season_start_day = 50#c(50, 160)  # days from start of sim
  )
  # Generate LHS samples
  A <- randomLHS(n_param_sets, 4)
  # Scale to parameter ranges
  params_df <- data.frame(
    max_SMC_kill_rate = qunif(A[,1], param_ranges$max_SMC_kill_rate[1], param_ranges$max_SMC_kill_rate[2]),
    lambda = qunif(A[,2], param_ranges$lambda[1], param_ranges$lambda[2]),
    kappa = qunif(A[,3], param_ranges$kappa[1], param_ranges$kappa[2])#,
    # season_start_day = round(qunif(A[,4], param_ranges$season_start_day[1], param_ranges$season_start_day[2]))
  )
  params_df$sim_id <- paste0('parameter_set_', rownames(params_df),"_", country_to_run, "_", sim_allow_superinfections)
  
  prob_bite_generic <- readRDS(paste0(path, '/archive/fit_rainfall/20251009-144330-1d355186/prob_bite_generic.rds'))
  prob_bite_generic$prob_infectious_bite = 0.15
  p_bitevector <- calc_lagged_vectors(prob_bite_generic, 0, burnints = burnints) # no lagged values
  
  params_df$lag_p_bite <- 0
  params_df$season_start_day <- 50
  # SMC delivery
  params_df <- params_df %>%
    rowwise() %>%
    mutate(smc_dose_days = list(c(seq(season_start_day, season_start_day + 120 - 1, 30),
                                  seq(season_start_day + 365, season_start_day + 365 + 120 - 1, 30),
                                  seq(season_start_day + 365*2, season_start_day  + 365*2 + 120 - 1, 30)))) %>%
    ungroup()
  
  parameters_df <- params_df %>%
    mutate(
      p_bite = purrr::map(lag_p_bite, ~p_bitevector[[paste0("lag_", .x)]])
    )
  
  # Make list of parameters instead of df
  params_list <- split(parameters_df, seq(nrow(parameters_df))) #%>%
    # lapply(function(df_row) {
    #   l <- as.list(df_row)
    #   # flatten any one-element list-columns if you want
    #   l$smc_dose_days <- df_row$smc_dose_days[[1]]
    #   l$p_bite <- df_row$p_bite[[1]]
    #   l
    # })#%>% lapply(as.list)
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
      
      TRUE
    })
    
    parallel::clusterExport(cl, c("params_list", "metadata_df", "base_inputs", "gen_bs",
                                  "n_particles", "n_threads", "burnints", "threshold", "tstep",
                                  "t_liverstage", "country_to_run", "VB", "divide"),
                            envir = environment())
    
    results2 <- parallel::clusterApply(cl,
                                       params_list,
                                       function(params_row) {
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
                                         eff_daily <- calc_smc_efficacy(o$infection_records, 
                                                                        params_row, 
                                                                        by_week = FALSE)
                                         eff$sim_id <- params_row$sim_id
                                         eff_daily$sim_id <- params_row$sim_id
                                         
                                         return(list(efficacy_weekly = eff,
                                                     efficacy_daily = eff_daily,
                                                     params = params_row))
                                       }
    )
    parallel::stopCluster(cl)
 }
 
  
  efficacy_weekly <- purrr::map_df(results2, 'efficacy_weekly')
  efficacy_daily <- purrr::map_df(results2, 'efficacy_daily')
  params <- purrr::map_df(results2, 'params') %>%
    select(-smc_dose_days, -p_bite)
 
  saveRDS(params, paste0(path, '/src/fit_smc/outputs/parameters_', Sys.Date(), '.rds'))
  saveRDS(efficacy_weekly, paste0(path, '/src/fit_smc/outputs/efficacy_weekly_', Sys.Date(), '.rds'))
  saveRDS(efficacy_daily, paste0(path, '/src/fit_smc/outputs/efficacy_daily_', Sys.Date(), '.rds'))
}


# Observed efficacy extracted from figure in SI of Hayley's paper 10.1016/S2214-109X(22)00416-8
# observed_efficacy <- read.csv(paste0(path, '/shared/smc_fits_hayley.csv'))
# 
# # Look at efficacy output 
# ggplot(eff %>% filter(weeks_since_smc < 10)) + 
#   geom_point(aes(x = weeks_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.1) + 
#   geom_line(aes(x = weeks_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.2) + 
#   ylim(c(-0.5, 1)) +
#   geom_line(data = observed_efficacy, aes(x = day_since_smc/7, y = efficacy)) +
#   geom_point(data = observed_efficacy, aes(x = day_since_smc/7, y = efficacy)) +
#   scale_x_continuous(breaks = seq(0,9,1)) +
#   theme_minimal() +
#   theme(legend.position = 'none')  
# 
# # Assume observed data has Normal measurement error
# calculate_likelihood <- function(observed, predicted, sigma = 0.1) {
#   # log-likelihood
#   sum(dnorm(observed, mean = predicted, sd = sigma, log = TRUE), na.rm = TRUE)
# }
# 
# likelihoods <- eff %>%
#   ungroup() %>%
#   mutate(day_since_smc = round(weeks_since_smc * 7,0)) %>%
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
#            weeks_since_smc < 10) %>%
#   ggplot(aes(x = weeks_since_smc*7, y = efficacy, group = sim_id, color = as.factor(sim_id))) +
#   geom_line(alpha = 0.5) +
#   geom_point(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), 
#              color = "black", inherit.aes = FALSE) +
#   geom_line(data = observed_efficacy, aes(x = day_since_smc, y = efficacy),
#             color = "black", inherit.aes = FALSE) +
#   labs(title = "Top 10 runs vs observed SMC efficacy",
#        y = "SMC efficacy", x = "Days since SMC") +
#   theme_minimal() + ylim(c(-0.5, 1))
# 
# # parameter_set_293_generic_TRUE is closest
# # max_SMC_kill_rate lambda kappa season_start_day sim_id                         lag_p_bite
# # <dbl>             <dbl>  <dbl> <dbl>            <chr>                          <dbl>
# # 6.41              44.2   0.425 150              parameter_set_293_generic_TRUE 0
# 
# top 10 
# max_SMC_kill_rate lambda kappa season_start_day sim_id                         lag_p_bite
# <dbl>  <dbl> <dbl>            <dbl> <chr>                               <dbl>
#   1             12.7   41.2  8.20                98 parameter_set_114_generic_TRUE          0
# 2              8.80  36.8  9.55                51 parameter_set_133_generic_TRUE          0
# 3              7.99   9.99 0.235               41 parameter_set_248_generic_TRUE          0
# 4              6.41  44.2  0.425              150 parameter_set_293_generic_TRUE          0
# 5             24.5   30.0  3.71                99 parameter_set_420_generic_TRUE          0
# 6             20.0   36.1  7.71                89 parameter_set_545_generic_TRUE          0
# 7             21.0   32.1  4.15                37 parameter_set_802_generic_TRUE          0
# 8             16.0   33.5  8.39                74 parameter_set_867_generic_TRUE          0
# 9             14.4   33.5  9.87                94 parameter_set_934_generic_TRUE          0
# 10             22.7   31.0  9.05               123 parameter_set_993_generic_TRUE          0

# eff %>%
#   filter(sim_id %in% top_runs$sim_id &
#            weeks_since_smc < 10) %>%
#   ggplot() +
#   geom_line(aes(x = weeks_since_smc, y = inci_none, group = sim_id, color = as.factor(sim_id)), 
#             alpha = 0.8, linetype = 2) +
#   geom_line(aes(x = weeks_since_smc, y = inci_smc, group = sim_id, color = as.factor(sim_id)), 
#             alpha = 0.8) +
#   labs(y = "incidence", x = "weeks since SMC") +
#   theme_minimal() 
