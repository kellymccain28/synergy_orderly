# Script to do a generic cohort simulation 
# almost the same as sim_cohort_grid but just separated to make clear 
sim_cohort_generic <- function(trial_ts = 365*3, 
                               treatment_prob = 0.9, # default is 1 (which gives children prophylaxis)
                               season_start_day = 137, # default is 137 to start on August 15 (days since April 1)
                               vax_day = 75, # default is mid june for 3rd dose 
                               threshold = 5000, # default is 5000 parasites per microL
                               N = 2000,
                               country_to_run = 'generic',
                               season = 'seasonal',
                               n_param_sets,
                               path = "R:/Kelly/synergy_orderly/",
                               notes,# notes are to write down the specifics of the runs 
                               get_parasit = FALSE){ 
  library(odin2)
  library(ggplot2)
  library(dust2)
  library(tidyverse)
  library(lhs)
  library(cyphr)
  
  source('R:/Kelly/synergy_orderly/shared/cohort_sim_utils.R')
  source('R:/Kelly/synergy_orderly/shared/helper_functions.R')
  source("R:/Kelly/synergy_orderly/shared/rtss.R")
  source("R:/Kelly/synergy_orderly/shared/likelihood.R")
  source("R:/Kelly/synergy_orderly/src/fit_smc/calculate_efficacy_likelihood.R")

  # sim_pars <- paste0("_treat_", treatment_prob, "start_", season_start_day, "threshold", threshold)

    # Load the within-host model 
  gen_bs <- odin2::odin("R:/Kelly/synergy_orderly/shared/smc_rtss.R")
  message('got model')

  # set base parameters 
  n_particles = 1L
  n_threads = 1L
  burnints = 50#90
  # threshold = 5000
  tstep = 1
  t_liverstage = 7 # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC267587/
  country_short = str_sub(country_to_run, 1, 1)
  VB = 1e6
  divide = if(tstep == 1) 2 else 1
  # Day of intervention (0 = start of follow-up; - values are before follow-up; + values after follow-up) - where 0 is also end of burnin
  # these get converted later to the correct directon - i.e. vaccine before follow-up will be +, smc before follow up will be -
  # vax_day is the 3rd primary dose (when we assume that efficacy begins)
  # vax_day = 75 # unlike the model sim, this is in days (not timesteps), ~75 days is mid-June which in the generic sim is just as the season is starting 
  # N = 4000
  
  treatment_probability = treatment_prob # in trial, everyone who was diagnosed with clincial malaria was treated 
  successful_treatment_probability = 0.9 # AL treatment protects up to 90% for 12 days SI of Commun. 5:5606 doi: 10.1038/ncomms6606 (this is what is in malsim)
  
  # Set up base inputs (these don't vary across parameter sweep)
  base_inputs <- list(
    trial_timesteps = trial_ts,#365*3,
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
  
  # Create parameter grid
  set.seed(123)
  
  params_df <- params_df <- data.frame(
    max_SMC_kill_rate = rep(2.37, n_param_sets),
    lambda = rep(18.5, n_param_sets),
    kappa = rep(0.337, n_param_sets),
    alpha_ab = rep(1.74, n_param_sets),# from 01-23: rep(1.66, n_param_sets), #rep(1.77, n_param_sets),#rep(1.38, n_param_sets),
    beta_ab = rep(4.69, n_param_sets),# from 01-23: rep(3.45, n_param_sets),  #rep(2.63, n_param_sets),#rep(5.83, n_param_sets)
    vmin = rep(0.00259, n_param_sets)# from 01-23: rep(0.00311, n_param_sets) #rep(0.000513, n_param_sets)
  )
  params_df$sim_id <- paste0('parameter_set_', rownames(params_df),"_", country_to_run, "_", treatment_probability)
  
  # probability of a bite is used in the cohort simulation and so is in 1-day timesteps
  # prob_bite_generic <- readRDS(paste0(path, 'archive/fit_rainfall/20251009-144330-1d355186/prob_bite_generic.rds'))
  if(season == 'seasonal'){
    prob_bite_generic <- readRDS("R:/Kelly/synergy_orderly/archive/fit_rainfall/20251209-213103-7f252a05/prob_bite_generic.rds")#with EIR of 30 and seasonal instead of seasonal(03/12/25)
  } else if (season == 'perennial'){
    prob_bite_generic <- readRDS("R:/Kelly/synergy_orderly/archive/fit_rainfall/20251209-213103-7f252a05/prob_bite_generic_perennial.rds")
  } else if (season == 'constant') {
    prob_bite_generic <- readRDS("R:/Kelly/synergy_orderly/archive/fit_rainfall/20251209-213103-7f252a05/prob_bite_generic_perennial.rds")
    prob_bite_generic$prob_infectious_bite <- 0.15 #constant
  }
  message('stopped after getting prob bite')
  p_bitevector <- calc_lagged_vectors(prob_bite_generic, 0, burnints = burnints) # no lagged values 
  
  params_df$lag_p_bite <- 0
  
  # SMC delivery 
  # season_start_day <- 137 # days between April 1 to August 15, so if the sim starts on 1 April, then SMC is delivered for 4 months Aug,Sep,Oct,Nov
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
    # could try to reduce effective coverage of SMC -- in the fitting
    # mutate(SMC = ifelse(arm == 'smc' & runif(N) > 0.9, 0, SMC)) %>%
    mutate(t_to_boost1 = 365,
           t_to_boost2 = 730,
           country = country_to_run) %>%
    mutate(rid_original =paste0(country_short, sprintf("%04d", rid)),
           country = 'generic',
           v1_date = as.Date('2017-04-01') + vax_day - 60) %>%# before was just april1 (start of), but now the first dose is dependent on the 3rd(specified by vax_day), 3rd dose is ~60 days after first
    mutate(fu_end_date = as.Date('2020-03-31'),
           fu_end_day = as.Date('2020-03-31') - as.Date('2017-04-01'))
  
  output_dir = paste0(path, 'src/sim_cohort_generic/outputs/outputs_', Sys.Date())
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
  # Save parameter grid
  saveRDS(parameters_df, paste0(output_dir, "/parameter_grid.rds"))
  saveRDS(base_inputs, paste0(output_dir, "/base_inputs.rds"))
  saveRDS(metadata_df, paste0(output_dir, "/metadata_df.rds"))
  writeLines(notes, paste0(output_dir, "/sim_notes.txt"))
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
                                                    # allow_superinfections = TRUE,
                                                    return_parasitemia = get_parasit,
                                                    save_outputs = FALSE)
                         message('finished simulation')
                         o$infection_records$sim_id <- params_row$sim_id
                         
                         return(list(infection_records = o$infection_records, 
                                     # parasitemia = o$parasitemia_data,
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
      source("R:/Kelly/synergy_orderly/shared/format_model_output.R")
      source("R:/Kelly/synergy_orderly/shared/get_incidence.R")
      
      TRUE
    })
    
    parallel::clusterExport(cl, c("params_list", "metadata_df", "base_inputs", "gen_bs",
                                  "n_particles", "n_threads", "burnints", "threshold", "tstep",
                                  "t_liverstage", "country_to_run", "VB", "divide", "output_dir"),
                            envir = environment())
    
    results2 <- parallel::clusterApply(cl,
                                       params_list,
                                       function(params_row) {
                                         o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                                                                    metadata_df,
                                                                    base_inputs,
                                                                    output_dir = 'R:/Kelly/src/sim_cohort_generic/outputs',
                                                                    # allow_superinfections = TRUE,
                                                                    return_parasitemia = get_parasit,
                                                                    save_outputs = FALSE)
                                         message('finished simulation')
                                         o$infection_records$sim_id <- params_row$sim_id
                                         
                                         # remove any infections that occurred within 7 days 
                                         infs <- o$infection_records %>% 
                                           group_by(rid) %>%
                                           arrange(rid, detection_day) %>%
                                           mutate(previous_detday = lag(detection_day),
                                                  diff = detection_day - previous_detday) %>%
                                           filter(diff > 7 | is.na(diff)) %>% select(-diff, -previous_detday)
                                         
                                         infs_formatted <- format_model_output(model_data = infs, 
                                                                               cohort = 'generic', 
                                                                               start_cohort = as.Date('2017-04-01'),
                                                                               simulation = params_row$sim_id)
                                         
                                         inci <- get_incidence(df_children = metadata_df,
                                                       casedata = infs_formatted) %>%
                                           mutate(sim_id = params_row$sim_id)
                                         
                                         if(get_parasit){
                                           out <- list(infection_records_formatted = infs_formatted, 
                                              incidence_df = inci,
                                              parasitemia = o$parasitemia_data,
                                              params = params_row)
                                         } else {
                                           out <- list(infection_records_formatted = infs_formatted, 
                                                       incidence_df = inci,
                                                       params = params_row)
                                           }
                                         
                                         return(out)
                                       }
    )
    parallel::stopCluster(cl)
  }
  # Determine number of batches
  n_results <- length(results2)
  batch_size <- 8
  n_batches <- ceiling(n_results / batch_size)
  
  # Loop through batches and save each one
  for(batch in 1:n_batches) {
    # Calculate start and end indices for this batch
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, n_results)
    
    # Extract for this batch
    infectionsformatted_batch <- purrr::map_df(results2[start_idx:end_idx], "infection_records_formatted")#"infection_records")
    incidence_batch <- purrr::map_df(results2[start_idx:end_idx], "incidence_df")#"infection_records")
    
    # Save this batch
    saveRDS(infectionsformatted_batch, paste0(output_dir, '/infs_formatted_batch', batch, '.rds'))
    saveRDS(incidence_batch, paste0(output_dir, '/incidence_batch', batch, '.rds'))
  }
  
  if(get_parasit){
    parasitemia <- purrr::map_df(results2, 'parasitemia')
    saveRDS(parasitemia, paste0(output_dir, "/parasitemia.rds"))
  }
  
  # infectionrecords <- purrr::map_df(results2, "infection_records")
  params <- purrr::map_df(results2, 'params')
  # saveRDS(infectionrecords, paste0(output_dir, '/infection_records.rds'))
  saveRDS(params, paste0(output_dir, "/parameter_df.rds"))
  
}