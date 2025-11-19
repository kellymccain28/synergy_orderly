# Script to do a generic cohort simulation 
# almost the same as sim_cohort_grid but just separated to make clear 
sim_cohort_generic <- function(trial_ts = 365*3, 
                               sim_allow_superinfections = TRUE,
                               country_to_run = 'generic',
                               n_param_sets,
                               path = "R:/Kelly/synergy_orderly/"){
  library(odin2)
  library(ggplot2)
  library(dust2)
  library(tidyverse)
  library(lhs)
  library(cyphr)
  
  # orderly_strict_mode()
  # orderlyparams <- orderly_parameters(#N = NULL, # size of cohort population 
  #   trial_ts = NULL,# trial timesteps in cohort simulation (inte)
  #   sim_allow_superinfections = NULL, # TRUE or FALSE 
  #   country_to_run = NULL, # BF or Mali, or if generic, then 'generic' which means that metadata_df is different.
  #   n_param_sets = NULL)  
  source('R:/Kelly/synergy_orderly/shared/cohort_sim_utils.R')
  source('R:/Kelly/synergy_orderly/shared/helper_functions.R')
  source("R:/Kelly/synergy_orderly/shared/rtss.R")
  source("R:/Kelly/synergy_orderly/shared/likelihood.R")
  source("R:/Kelly/synergy_orderly/src/fit_smc/calculate_efficacy_likelihood.R")
  # orderly_shared_resource("rtss.R",
  #                         "smc_rtss.R",
  #                         "helper_functions.R",
  #                         "format_model_output.R",
  #                         "get_incidence.R",
  #                         "get_cox_efficacy.R",
  #                         "analyse_model_output.R",
  #                         "likelihood.R",
  #                         "cohort_sim_utils.R")
  # 
  # orderly_dependency(name = 'fit_rainfall',
  #                    "latest()",
  #                    c("prob_bite_BFA.rds",
  #                      "prob_bite_MLI.rds",
  #                      'prob_bite_generic.rds'
  #                    ))
  
  # # Source antibody function
  # source("rtss.R")
  # # Source helper functions
  # source("helper_functions.R")
  # Load the within-host model 
  gen_bs <- odin2::odin("R:/Kelly/synergy_orderly/shared/smc_rtss.R")
  message('got model')
  # Source the utils functions 
  # source("cohort_sim_utils.R")
  # # SOurce processing functions 
  # source("likelihood.R")
  # source('get_cox_efficacy.R')
  # source("format_model_output.R")
  # source("get_incidence.R")
  # source("analyse_model_output.R")
  
  # set base parameters 
  n_particles = 1L
  n_threads = 1L
  burnints = 50#90
  threshold = 5000
  tstep = 1
  t_liverstage = 8
  country_short = str_sub(country_to_run, 1, 1)
  VB = 1e6
  divide = if(tstep == 1) 2 else 1
  # Day of intervention (0 = start of follow-up; - values are before follow-up; + values after follow-up) - where 0 is also end of burnin
  # these get converted later to the correct directon - i.e. vaccine before follow-up will be +, smc before follow up will be -
  # vax_day is the 3rd primary dose (when we assume that efficacy begins)
  vax_day = -10 # unlike hte model sim, this is in days (not timesteps)
  N = 400
  
  # Set up base inputs (these don't vary across parameter sweep)
  base_inputs <- list(
    trial_timesteps = trial_ts,#365*3,
    burnin = burnints,
    threshold = threshold,
    VB = VB,
    tstep = tstep,
    t_liverstage = t_liverstage,
    country = country_to_run,
    country_short = country_short
  )
  
  # Create parameter grid
  set.seed(123)
  
  params_df <- params_df <- data.frame(
    max_SMC_kill_rate = rep(3, n_param_sets),
    lambda = rep(13.08, n_param_sets),
    kappa = rep(0.43, n_param_sets)
  )
  params_df$sim_id <- paste0('parameter_set_', rownames(params_df),"_", country_to_run, "_", sim_allow_superinfections)
  
  # probability of a bite is used in the cohort simulation and so is in 1-day timesteps
  prob_bite_generic <- readRDS(paste0(path, 'archive/fit_rainfall/20251009-144330-1d355186/prob_bite_generic.rds'))
  message('stoppedafter getting prob bite')
  p_bitevector <- calc_lagged_vectors(prob_bite_generic, 0, burnints = burnints) # no lagged values 
  
  params_df$lag_p_bite <- 0
  
  # SMC delivery 
  season_start_day <- 122 # days between April 1 to August 1, so if the sim starts on 
  season_length <- 120
  params_df <- params_df %>%
    rowwise() %>%
    mutate(smc_dose_days = list(c(seq(season_start_day, season_start_day + season_length - 1, 30),
                                  seq(season_start_day + 366, season_start_day + 366 + season_length - 1, 30),
                                  seq(season_start_day + 365*2, season_start_day  + 365*2 + season_length - 1, 30))))
  
  # Add to parameter dataframe by matching lag values (always 0 for generic)
  parameters_df <- params_df %>%
    mutate(
      p_bite = map(lag_p_bite, ~p_bitevector[[paste0("lag_", .x)]])
    )
  
  # Make list of parameters instead of df 
  params_list <- split(parameters_df, seq(nrow(parameters_df)))
  
  metadata_df <- data.frame(
    rid = 1:N,
    vaccination_day = vax_day,#sample(-1:0, N, replace = TRUE), the model takes abs(days before start of follow-up)
    PEV = rbinom(N, 1, 0.5),#rep(0,N),#
    SMC = c(rep(1, N/2), rep(0, N/2))#rbinom(N, 1, 0.5)
  ) %>%
    mutate(arm = case_when(
      PEV == 1 & SMC == 1 ~ 'both',
      PEV == 1 & SMC == 0 ~ 'rtss',
      PEV == 0 & SMC == 1 ~ 'smc',
      TRUE ~ 'none')) %>%
    mutate(t_to_boost1 = 365,
           t_to_boost2 = 730,
           country = country_to_run) %>%
    mutate(rid_original =paste0(country_short, sprintf("%04d", rid)),
           country = 'generic',
           v1_date = as.Date('2017-04-01'))
  
  output_dir = paste0(path, 'src/sim_cohort_generic/outputs/outputs_', Sys.Date())
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # Save parameter grid
  saveRDS(parameters_df, file.path(output_dir, "parameter_grid.rds"))
  saveRDS(base_inputs, file.path(output_dir, "base_inputs.rds"))
  saveRDS(metadata_df, file.path(output_dir, "metadata_df.rds"))
  message('stopped here after saving')
  
  # Run simulation
  cluster_cores <- Sys.getenv("CCP_NUMCPUS")
  # if (cluster_cores == "") {
  #   cluster_cores <- 8
  # }
  
  if (cluster_cores == "") {
    message("running in serial (on a laptop?)")
    message("Running ", nrow(params_df), " simulations sequentially")
    
    results2 <- lapply(params_list,
                       function(params_row){
                         o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                                                    metadata_df,
                                                    base_inputs,
                                                    output_dir = 'R:/Kelly/src/sim_cohort_generic/outputs',
                                                    allow_superinfections = TRUE,
                                                    return_parasitemia = FALSE,
                                                    save_outputs = FALSE)
                         message('finished simulation')
                         o$infection_records$sim_id <- params_row$sim_id
                         
                         return(list(infection_records = o$infection_records, 
                                     parasitemia = o$parasitemia_data,
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
                                                                    output_dir = 'R:/Kelly/src/sim_cohort_generic/outputs',
                                                                    allow_superinfections = TRUE,
                                                                    return_parasitemia = TRUE,
                                                                    save_outputs = FALSE)
                                         message('finished simulation')
                                         o$infection_records$sim_id <- params_row$sim_id
                                         
                                         return(list(infection_records = o$infection_records, 
                                                     parasitemia = o$parasitemia_data,
                                                     params = params_row))
                                       }
    )
    parallel::stopCluster(cl)
  }
  
  saveRDS(results2, paste0(output_dir, '/sim_results.rds'))
}
# sim_results <- readRDS("R:/Kelly/synergy_orderly/src/sim_cohort_generic/outputs/outputs_2025-11-19/sim_results.rds")
# 
# infectionrecords <- purrr::map_df(sim_results, "infection_records")
# params <- purrr::map_df(sim_results, 'params')
# metadata_df <- readRDS("R:/Kelly/synergy_orderly/src/sim_cohort_generic/outputs/outputs_2025-11-19/metadata_df.rds")
# 
# # saveRDS(params, paste0(output_dir, '/parameters_', Sys.Date(), '.rds'))
# # saveRDS(infectionrecords, paste0(output_dir, '/infectionrecords_', Sys.Date(), '.rds'))
# formattedinfrecords <- lapply(params$sim_id, function(x){
#   format_model_output(model_data = infectionrecords,
#                       cohort = 'generic',
#                       simulation = x)})
# all <- bind_rows(formattedinfrecords)
# testinci <- lapply(params$sim_id, function(x){
#   aa <- all %>% filter(sim_id == x)
# 
#   get_incidence(df_children = metadata_df,
#                           casedata = aa) %>%
#     mutate(sim_id = x)
# })
# 
# ggplot(inci %>% filter(year < 2018)) +
#   geom_line(aes(x = yearmonth, y = incidence_per_1000pm, color = arm)) +
#   # geom_vline(xintercept = zoo::as.yearmon(as.Date(unlist(all$smc_dose_days[1][1:4])-90, origin = '2017-01-01'))) +
#   facet_wrap(~sim_id)
# # ok here it was delivered on day 212 from jan 1, but if we remove 90 days, which is the difference between jan 1 and april 1