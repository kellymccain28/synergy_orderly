# Script to fit spline for prob of a bite in cohort of the trial simulation 
optimise_sim_trial_cohort <- function(trial_ts = 365*3, 
                                      treatment_prob = 1, # default is 1 (which gives children prophylaxis)
                                      threshold = 5000, # default is 5000 parasites per microL
                                      country_to_run, # should be BF or Mali
                                      n_param_sets,
                                      path = "R:/Kelly/synergy_orderly/",
                                      get_parasit = FALSE,
                                      notes){
  
  # almost the same as sim_cohort_generic and sim_trial_cohort but separate for clarity 
    library(odin2)
    library(ggplot2)
    # library(dust)
    library(tidyverse)
    library(lhs)
    library(zoo)
    
    setwd(path)
    output_dir = paste0(path, 'src/sim_trial_cohort/outputs_fitting/outputs_', Sys.Date())
    
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
    source("R:/Kelly/synergy_orderly/shared/format_model_output.R")
    source("R:/Kelly/synergy_orderly/shared/get_incidence.R")
    source("R:/Kelly/synergy_orderly/src/sim_trial_cohort/optimisation_helpers.R")
    
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
    
    params_df <- data.frame(#params_df %>%
      # mutate(
      max_SMC_kill_rate = rep(2.37, n_param_sets),
      lambda = rep(18.5, n_param_sets),
      kappa = rep(0.337, n_param_sets),
      alpha_ab = rep(1.74, n_param_sets),
      beta_ab = rep(4.69,  n_param_sets),
      vmin = rep(0.00259, n_param_sets))
    params_df$lag_p_bite <- 0
    
    params_df$sim_id <- paste0('parameter_set_', rownames(params_df),"_", country_to_run)
    
    # probability of a bite is used in the cohort simulation and so is in 1-day timesteps
    # if(country_to_run == 'BF') {
    #   prob_bite <- readRDS('R:/Kelly/synergy_orderly/archive/fit_rainfall/20251209-213103-7f252a05/prob_bite_BFA.rds') 
    #   prob_bite$prob_infectious_bite <- prob_bite$prob_infectious_bite 
    # } else if(country_to_run == 'Mali') {
    #   prob_bite <- readRDS('R:/Kelly/synergy_orderly/archive/fit_rainfall/20251209-213103-7f252a05/prob_bite_MLI.rds') 
    #   prob_bite$prob_infectious_bite <- prob_bite$prob_infectious_bite
    # }
    # p_bitevector <- calc_lagged_vectors(prob_bite, unique(params_df$lag_p_bite), burnints = burnints) # no lagged values 
    # 
    # # Add to parameter dataframe by matching lag values 
    # parameters_df <- params_df %>%
    #   mutate(
    #     p_bite = map(lag_p_bite, ~p_bitevector[[paste0("lag_", .x)]])
      # )
    
    # Load saved spline setup
    if(country_to_run =='BF'){
      setup <- readRDS(paste0(path, "src/sim_trial_cohort/spline_setupBF.rds"))
    } else if(country_to_run == 'Mali'){
      setup <- readRDS(paste0(path, "src/sim_trial_cohort/spline_setupMali.rds"))
    }
    
    
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
      mutate(rid_original = paste0(country_short, rid))
    
    # Join synthetic no intervention and real trial kids together
    metadata_df <- full_join(metadata_df, noint_arm)
    
    # Save parameter grid
    saveRDS(params_df, paste0(output_dir, "/parameter_grid.rds"))
    saveRDS(base_inputs, paste0(output_dir, "/base_inputs.rds"))
    saveRDS(metadata_df, paste0(output_dir, "/metadata_df.rds"))
    writeLines(notes, paste0(output_dir, "/sim_notes.txt"))
    message('stopped here after saving')
    
    # Get incidence to compare to depending on the country run
    if(country_to_run == 'BF'){
      incidence_trial <- readRDS('R:/Kelly/synergy_orderly/archive/trial_results/20260219-104643-cb128f65/monthly_incidence_trial_BF.rds')#readRDS("R:/Kelly/synergy_orderly/archive/trial_results/20260127-164922-12477275/monthly_incidence_trial.rds")
    } else if(country_to_run == 'Mali'){
      incidence_trial <- readRDS('R:/Kelly/synergy_orderly/archive/trial_results/20260219-104643-cb128f65/monthly_incidence_trial_Mali.rds')
    }
    
    # # To resume if crashed
    # if(file.exists("optim_checkpoint.rds")) {
    #   cat("Resuming from checkpoint...\n")
    #   optimization_state <- readRDS("optim_checkpoint.rds")
    #   # Use the last coefficients as starting point
    #   starting_coefs <- optimization_state$best_coefs
    # } else {
    #   # Fresh start
    #   starting_coefs <- your_initial_coefs
    # }
    
    # Run simulation
    cluster_cores <- Sys.getenv("CCP_NUMCPUS")
    
    if (cluster_cores == "") {
      # Serial version
      results2 <- lapply(1:nrow(params_df), function(i) {
        fit_prob_bite_spline(
          params_row = params_df[i,],
          metadata_df = metadata_df,
          base_inputs = base_inputs,
          output_dir = output_dir,
          incidence_trial = incidence_trial,
          setup = setup)
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
        source("R:/Kelly/synergy_orderly/src/sim_trial_cohort/optimisation_helpers.R")
        
        TRUE
      })
      
      parallel::clusterExport(cl, c("params_df", "metadata_df", "base_inputs", "gen_bs",
                                    "n_particles", "n_threads", "burnints", "threshold", "tstep",
                                    "t_liverstage", "country_to_run", "VB", "divide", "output_dir",
                                    "incidence_trial", "setup"),
                              envir = environment())
      
      # Run in parallel (one per parameter set)
      results2 <- parallel::parLapply(cl, 1:nrow(params_df), function(i) {
        fit_prob_bite_spline(
          params_row = params_df[i,],
          metadata_df = metadata_df,
          base_inputs = base_inputs,
          output_dir = output_dir,
          incidence_trial = incidence_trial,
          setup = setup,
          country_to_run = country_to_run
        )
      })                               
                      
      saveRDS(results2, 'final_results.rds')               
      parallel::stopCluster(cl)
      
      
    }
  }
  
  