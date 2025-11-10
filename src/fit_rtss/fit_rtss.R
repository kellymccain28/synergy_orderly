run_fit_rtss <- function(path = "R:/Kelly/synergy_orderly",
                        n_param_sets,
                        N = 1200){
  path = "R:/Kelly/synergy_orderly"
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
  source(paste0(path, "/src/sim_cohort_grid/cohort_sim_utils.R"))
  # SOurce processing functions
  source(paste0(path, "/shared/likelihood.R"))
  
  trial_ts = 365# trial timesteps in cohort simulation (inte)
  sim_allow_superinfections = TRUE # TRUE or FALSE
  country_to_run = 'generic - vaccine'
  country_short = 'g'
  n_param_sets = n_param_sets
  N = N
  vax_day = -1 # unlike the model sim, this is in days (not timesteps)
  
  n_particles = 1L
  n_threads = 1L
  burnints = 70
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
  
  params_df <- data.frame(
    max_SMC_kill_rate = rep(0, n_param_sets),
    lambda = rep(0, n_param_sets),
    kappa = rep(0, n_param_sets)
  )
  params_df$sim_id <- paste0('parameter_set_', rownames(params_df),"_", country_to_run, "_", sim_allow_superinfections)
  
  prob_bite_generic <- readRDS(paste0(path, '/archive/fit_rainfall/20251009-144330-1d355186/prob_bite_generic.rds'))
  prob_bite_generic$prob_infectious_bite = 0.3
  p_bitevector <- calc_lagged_vectors(prob_bite_generic, 0, burnints = burnints) # no lagged values
  
  params_df$lag_p_bite <- 0
  params_df$season_start_day <- 50
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
           country = 'generic - vaccine',
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
                                                    output_dir = 'R:/Kelly/src/fit_rtss/outputs',
                                                    allow_superinfections = TRUE,
                                                    return_parasitemia = FALSE,
                                                    save_outputs = FALSE)
                         message('finished simulation')
                         o$infection_records$sim_id <- params_row$sim_id
                         eff <- calc_rtss_efficacy(o$infection_records)# eff <- calc_smc_efficacy(o$infection_records,
                         #                          params_row, 
                         #                          by_week = TRUE)
                         # eff_daily <- calc_smc_efficacy(o$infection_records, 
                         #                                params_row, 
                         #                                by_week = FALSE)
                         # eff$sim_id <- params_row$sim_id
                         # eff_daily$sim_id <- params_row$sim_id
                         
                         # return(list(efficacy_weekly = eff,
                         #             efficacy_daily = eff_daily,
                         #             params = params_row))
                         return(list(infection_records = o$infection_records, 
                                     efficacy_weekly = eff,
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
                                                                    output_dir = 'R:/Kelly/src/fit_rtss/outputs',
                                                                    allow_superinfections = TRUE,
                                                                    return_parasitemia = FALSE,
                                                                    save_outputs = FALSE)
                                         message('finished simulation')
                                         o$infection_records$sim_id <- params_row$sim_id
                                         
                                         eff <- calc_rtss_efficacy(o$infection_records)
                                         # eff <- calc_smc_efficacy(o$infection_records,
                                         #                          params_row, 
                                         #                          by_week = TRUE)
                                         # eff_daily <- calc_smc_efficacy(o$infection_records, 
                                         #                                params_row, 
                                         #                                by_week = FALSE)
                                         # eff$sim_id <- params_row$sim_id
                                         # eff_daily$sim_id <- params_row$sim_id
                                         
                                         # return(list(efficacy_weekly = eff,
                                         #             efficacy_daily = eff_daily,
                                         #             params = params_row))
                                         return(list(infection_records = o$infection_records, 
                                                     efficacy_weekly = eff,
                                                     params = params_row))
                                       }
    )
    parallel::stopCluster(cl)
  }
  
  infectionrecords <- purrr::map_df(results2, "infection_records")
  efficacy <- purrr::map_df(results2, 'efficacy_weekly')
  params <- purrr::map_df(results2, 'parameters')
  
  saveRDS(params, paste0(path, '/src/fit_rtss/outputs/parameters_', Sys.Date(), '.rds'))
  saveRDS(infectionrecords, paste0(path, '/src/fit_rtss/outputs/infectionrecords_rtss_', Sys.Date(), '.rds'))
  saveRDS(efficacy, paste0(path, '/src/fit_rtss/outputs/efficacy_rtss_', Sys.Date(), '.rds'))
}


# Efficacy against infection from White 2015
# ts <- seq(1,365*3)
# phases <- ifelse(ts < 366, 1,
#                  ifelse(ts < 729, 2, 3))
# csp <- list()
# for(i in 1:50){
#   csp[[i]] <- antibody_titre(t = ts,
#                      phase = phases,
#                      peak1 = c(621,0.35) ,
#                      peak2 = c(277, 0.35),
#                      peak3 = c(277,0.35),
#                      duration1 = c(45,16),
#                      duration2 = c(591,245),
#                      rho1 = c(2.37832, 1.00813),
#                      rho2 = c(1.034, 1.027),
#                      rho3 = c(1.034, 1.027))
# }
# csp_df <- data.frame(do.call(rbind, csp))
# colnames(csp_df) <- ts
# csp_df$sim <- 1:50
# 
# csp_long <- pivot_longer(csp_df, cols = -sim,
#                          names_to = "time", values_to = "titre")
# csp_long$time <- as.numeric(csp_long$time)
# 
# # Plot
# ggplot(csp_long, aes(x = time, y = titre, group = sim)) +
#   geom_line(alpha = 0.3) +
#   theme_minimal() + scale_y_log10() +
#   labs(x = "Time", y = "Antibody Titre")
# 
# ve <- lapply(csp, vaccine_eff)
# ve_df <- data.frame(do.call(rbind, ve))
# colnames(ve_df) <- ts
# ve_df$sim <- 1:50
# 
# ve_long <- pivot_longer(ve_df, cols = -sim,
#                          names_to = "time", values_to = "ve")
# ve_long$time <- as.numeric(ve_long$time)
# ve_long <- ve_long %>% ungroup() %>%
#   mutate(weeks_since_rtss = floor(time/7)) %>%
#   filter(weeks_since_rtss < 52) %>%
#   group_by(weeks_since_rtss, sim) %>%
#   summarize(ve_inf = mean(ve))
# 
# # Plot
# ggplot(ve_long, aes(x = weeks_since_rtss, y = ve_inf, group = sim)) +
#   geom_line(alpha = 0.3) +
#   theme_minimal() +
#   labs(x = "Weeks since RTSS", y = "Efficacy")

# ve_inf <- vaccine_eff(csp)
# plot(ve_inf[0:52])
# ve_inf <- ve_inf[seq_along(ve_inf) %% 7 == 0]
# ve_inf_df <- data.frame(#weeks_since_rtss = seq(0,66),
#                         days_since_rtss = ts,
#                         ve_inf = ve_inf) %>%
#   mutate(weeks_since_rtss = floor(days_since_rtss/7)) %>%
#   group_by(weeks_since_rtss) %>%
#   summarize(ve_inf = mean(ve_inf))
# plot(ve_inf_df)
# infectionrecords <- readRDS("R:/Kelly/synergy_orderly/src/fit_rtss/outputs/infectionrecords_rtss_2025-10-16.rds")
# infectionrecords <- readRDS("R:/Kelly/synergy_orderly/src/fit_rtss/outputs/infectionrecords_rtss_2025-10-17.rds")
# df <- infectionrecords
#
# eff <- calc_rtss_efficacy(df) %>%
#   left_join(ve_inf_df)
#
# eff1020 <- readRDS("R:/Kelly/synergy_orderly/src/fit_rtss/outputs/efficacy_rtss_2025-10-23.rds")
# infectionrecords_rtss1023 <- readRDS("R:/Kelly/synergy_orderly/src/fit_rtss/outputs/infectionrecords_rtss_2025-10-23.rds")
# eff <- eff1020 %>%#eff1017 %>%
#   left_join(ve_long)
# 
# eff %>% #filter(weeks_since_rtss< 365)  %>%
#   ggplot() +
#   geom_line(aes(x = weeks_since_rtss+8, y = ve_inf, group = sim),
#             alpha = 0.2, color = 'orange', linewidth = 1) +
#   geom_line(aes(x = weeks_since_rtss, y = efficacy, group = sim_id), color = '#962150', alpha = 0.2) +
#   ylim(c(0,1)) + xlim(c(0,60))+
#   theme_bw()
# ggsave(filename = 'outputs/efficacy_comparison_20251024_0.3_0_plus10.png')
