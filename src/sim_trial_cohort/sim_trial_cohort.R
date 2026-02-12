# Script to simulate the cohort of the trial simulation 

# almost the same as sim_cohort_generic but separate for clarity 
# major difference is inclusion of lag and scaler for seasonality 
sim_trial_cohort <- function(trial_ts = 365*3, 
                             treatment_prob = 0.9, # default is 1 (which gives children prophylaxis)
                             threshold = 5000, # default is 5000 parasites per microL
                             country_to_run, # should be BF or Mali
                             n_param_sets,
                             get_parasit = FALSE,
                             path = "R:/Kelly/synergy_orderly/",
                             notes){ # notes are to write down the specifics of the runs 
  library(odin2)
  library(ggplot2)
  # library(dust)
  library(tidyverse)
  library(lhs)
  setwd(path)
  output_dir = paste0(path, 'src/sim_trial_cohort/outputs/outputs_', Sys.Date())
  
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
  
  source('R:/Kelly/synergy_orderly/shared/cohort_sim_utils.R')
  source('R:/Kelly/synergy_orderly/shared/helper_functions.R')
  source("R:/Kelly/synergy_orderly/shared/rtss.R")
  source("R:/Kelly/synergy_orderly/shared/likelihood.R")
  source("R:/Kelly/synergy_orderly/src/sim_trial_cohort/compare_incidence.R")
  
  # Load the within-host model 
  gen_bs <- odin2::odin("R:/Kelly/synergy_orderly/shared/smc_rtss.R")
  message('loaded model')
  
  # Load data 
  # key <- cyphr::data_key()
  children <- readRDS('R:/Kelly/synergy_orderly/shared/children_deanonymised.rds')#cyphr::decrypt(readRDS("R:/Kelly/synergy_orderly/archive/clean_trial_data/20260127-163716-f220f556/data/children.rds"), key)
  
  # Get fu end date in days since april 1 2017 (beginning of trial)
  children$fu_end_day = round(children$fu_end_date - as.Date('2017-04-01'), 0)
  
  # set base parameters 
  n_particles = 1L
  n_threads = 1L
  burnints = 50#90
  tstep = 1
  t_liverstage = 7 # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC267587/
  country_short = str_sub(country_to_run, 1, 1)
  VB = 1e6
  divide = if(tstep == 1) 2 else 1
  
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
  A <- randomLHS(n = n_param_sets, k = 2) # n different sets of parameters, with k parameters to change
  A[,1] <- round(qunif(A[,1], -50, 50), 0) # lag in days of p of an infectious bite 
  A[,2] <- qunif(A[,2], 0.1, 1) # this is the scaling factor for the probability of infectious bite 
  colnames(A) <- c('lag_p_bite', 'p_bite_scaler')
  params_df <- as.data.frame(A)
  
  params_df <- params_df %>%
    mutate(
      max_SMC_kill_rate = 2.37,
      lambda = 18.5,
      kappa = 0.337,
      alpha_ab = 1.74,
      beta_ab = 4.69, 
      vmin = 0.00259)
    
  params_df$sim_id <- paste0('parameter_set_', rownames(params_df),"_", country_to_run)
  
  # probability of a bite is used in the cohort simulation and so is in 1-day timesteps
  if(country_to_run == 'BF') {
    prob_bite <- readRDS('R:/Kelly/synergy_orderly/archive/fit_rainfall/20251209-213103-7f252a05/prob_bite_BFA.rds') 
    prob_bite$prob_infectious_bite <- prob_bite$prob_infectious_bite 
  } else if(country_to_run == 'Mali') {
    prob_bite <- readRDS('R:/Kelly/synergy_orderly/archive/fit_rainfall/20251209-213103-7f252a05/prob_bite_MLI.rds') 
    prob_bite$prob_infectious_bite <- prob_bite$prob_infectious_bite
  }
  p_bitevector <- calc_lagged_vectors(prob_bite, unique(params_df$lag_p_bite), burnints = burnints) # no lagged values 
  
  # Add to parameter dataframe by matching lag values 
  parameters_df <- params_df %>%
    mutate(
      p_bite = map2(lag_p_bite, p_bite_scaler, 
                    ~ p_bitevector[[paste0("lag_", .x)]]$prob_lagged * .y)
    )
  
  # Make list of parameters instead of df 
  params_list <- split(parameters_df, seq(nrow(parameters_df)))
  
  # Make metadata data frame 
  # Create a metadata child data frame from trial data 
  metadata_df <- children %>%
    # filter to country we are running 
    filter(country == country_to_run) %>%
    mutate(start_of_fu = as.Date('2017-04-01'),  # this is the start of the cohort simulation. 
           start_to_v3 = v3_date - start_of_fu, # 
           fu_end_day = as.integer(fu_end_day),
           # vaccination day is the day of vaccination relative to the start of follow-up of cohort which is april 1, 2017
           # so a + value means that the 3rd dose of the vaccine was given x days after april 1, 2017
           # then we use this to find the number of days since vaccination within the cohrot sim loop to determine when efficacy begins
           vaccination_day = as.numeric(start_to_v3),
           PEV = ifelse(arm == 'smc', 0, 1),
           SMC = ifelse(arm == 'rtss', 0, 1),
           t_to_boost1 = as.numeric(boost1_date - v3_date),
           t_to_boost2 = as.numeric(boost2_date - v3_date)) %>%
    # Calculate timings of smc for each year - first, need to calculate days since april 1, 2017 (analogous to vaccination_day above)
    rowwise() %>%
    mutate(smc_dates = list(na.omit(c_across(ends_with('date_received')))),
           smc_dose_days = list(as.integer(smc_dates - start_of_fu)),
           smc_dose_days = if_else(arm == 'rtss', list(trial_ts+burnints), list(smc_dose_days))
    ) %>%
    # For missing vaccination days, make it very late 
    mutate(vaccination_day = ifelse(is.na(vaccination_day), 1350, vaccination_day)) %>%
    # for boosters, say if it is missing then second booster is much later so that it is after the follow-up time is over
    mutate(t_to_boost1 = ifelse(is.na(t_to_boost1), 1400, t_to_boost1),
           t_to_boost2 = ifelse(is.na(t_to_boost2), 1500, t_to_boost2)) %>%
    # update RIDs to be numeric 
    mutate(rid_original = rid) %>%
    group_by(country) %>%
    # make numeric id to use for darwing of numbers 
    mutate(rid = as.numeric(str_extract(rid, "[0-9]+"))) %>%
    arrange(rid) %>%
    mutate(rid = row_number())
  
  noint_arm <- data.frame(
    rid = seq(max(metadata_df$rid) + 1, max(metadata_df$rid) + 915),
    arm = 'none',
    sex = NA, 
    country = country_to_run, 
    start_of_fu = as.Date('2017-04-01'),
    start_to_v3 = NA, 
    vaccination_day = 1350, 
    PEV = 0, 
    SMC = 0, 
    t_to_boost1 = 1400, 
    t_to_boost2 = 1500, 
    smc_dates = NA, 
    smc_dose_days = NA,
    v1_date = as.Date('2017-04-01'), # this is because v1_date is used as start of FU
    fu_end_date = as.Date('2020-03-31'),
    fu_end_day = 1095,
    nsmc_received = 0
  ) %>%
    mutate(rid_original = paste0('B', rid))
  
  # Join synthetic no intervention and real trial kids together
  metadata_df <- full_join(metadata_df, noint_arm)
  
  # Save parameter grid
  saveRDS(parameters_df, paste0(output_dir, "/parameter_grid.rds"))
  saveRDS(base_inputs, paste0(output_dir, "/base_inputs.rds"))
  saveRDS(metadata_df, paste0(output_dir, "/metadata_df.rds"))
  writeLines(notes, paste0(output_dir, "/sim_notes.txt"))
  message('stopped here after saving')
  
  incidence_trial <- readRDS("R:/Kelly/synergy_orderly/archive/trial_results/20260127-164922-12477275/monthly_incidence_trial.rds")
  
  # Run simulation
  cluster_cores <- Sys.getenv("CCP_NUMCPUS")
  
  if (cluster_cores == "") {
    message("running in serial (on a laptop?)")
    message("Running ", nrow(params_df), " simulations sequentially")

    results2 <- lapply(params_list,
                       function(params_row){
                         o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                                                    metadata_df,
                                                    base_inputs,
                                                    output_dir = 'R:/Kelly/src/sim_trial_cohort/outputs',
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
                                                               cohort = country_to_run, 
                                                               start_cohort = as.Date('2017-04-01'),
                                                               simulation = params_row$sim_id)
                         
                         inci <- get_incidence(df_children = metadata_df,
                                               casedata = infs_formatted) %>%
                           mutate(sim_id = params_row$sim_id)
                         
                         metrics <- compare_incidence(incidence_model = inci, 
                                                      incidence_trial = incidence_trial)
                         
                         return(list(infection_records_formatted = infs_formatted, 
                                     incidence_df = inci,
                                     params = params_row,
                                     metrics = metrics))
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
      source("R:/Kelly/synergy_orderly/src/sim_trial_cohort/compare_incidence.R")
      
      TRUE
    })
    
    parallel::clusterExport(cl, c("params_list", "metadata_df", "base_inputs", "gen_bs",
                                  "n_particles", "n_threads", "burnints", "threshold", "tstep",
                                  "t_liverstage", "country_to_run", "VB", "divide", "output_dir",
                                  "incidence_trial"),
                            envir = environment())
    
    results2 <- parallel::clusterApply(cl,
                                       params_list,
                                       function(params_row) {
                                         o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                                                                    metadata_df,
                                                                    base_inputs,
                                                                    output_dir = 'R:/Kelly/src/sim_trial_cohort/outputs',
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
                                                                               cohort = country_to_run, 
                                                                               start_cohort = as.Date('2017-04-01'),
                                                                               simulation = params_row$sim_id)
                                         
                                         inci <- get_incidence(df_children = metadata_df,
                                                               casedata = infs_formatted) %>%
                                           mutate(sim_id = params_row$sim_id)
                                         
                                         metrics <- compare_incidence(incidence_model = inci, 
                                                                     incidence_trial = incidence_trial)
                                         
                                         return(list(infection_records_formatted = infs_formatted, 
                                                     incidence_df = inci,
                                                     params = params_row,
                                                     metrics = metrics))
                                       }
    )
    parallel::stopCluster(cl)
  
  
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
  
  metrics <- purrr::map_df(results2, 'metrics')
  saveRDS(metrics, paste0(output_dir, "/metrics.rds"))
  
  message('saved outputs at ', output_dir)
    }
}
