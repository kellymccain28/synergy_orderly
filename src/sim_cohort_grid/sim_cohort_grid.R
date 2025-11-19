library(odin2)
library(ggplot2)
library(dust2)
library(tidyverse)
library(mgcv)
library(umbrella)
library(lhs)
library(orderly2)
library(cyphr)
library(survival)
library(broom)
library(survminer)

key <- cyphr::data_key()

orderly_strict_mode()
orderlyparams <- orderly_parameters(#N = NULL, # size of cohort population 
  trial_ts = NULL,# trial timesteps in cohort simulation (inte)
  sim_allow_superinfections = NULL, # TRUE or FALSE 
  country_to_run = NULL, # BF or Mali, or if generic, then 'generic' which means that metadata_df is different.
  n_param_sets = NULL)  

orderly_shared_resource("rtss.R",
                        "smc_rtss.R",
                        "helper_functions.R",
                        "format_model_output.R",
                        "get_incidence.R",
                        "get_cox_efficacy.R",
                        "analyse_model_output.R",
                        "cohort_sim_utils.R",
                        "likelihood.R")

orderly_resource(c('run_process.R'))

orderly_dependency(name = 'fit_rainfall',
                   "latest()",
                   c("prob_bite_BFA.rds",
                     "prob_bite_MLI.rds",
                     'prob_bite_generic.rds'
                   ))

orderly_dependency(name = 'clean_trial_data',
                   "latest()",
                   c('data/children.rds'))

children <- cyphr::decrypt(readRDS("data/children.rds"), key)

orderly_dependency(name = 'trial_results',
                   "latest()",
                   c('monthly_incidence_trial.rds'))

# key <- cyphr::data_key()
monthly_inci_trial <- readRDS("monthly_incidence_trial.rds")

# Source antibody function
source("rtss.R")
# Source helper functions
source("helper_functions.R")
# Load the within-host model 
gen_bs <- odin2::odin("smc_rtss.R")
# Source the utils functions 
source("cohort_sim_utils.R")
# SOurce the function to run then process the runs
source('run_process.R')
# SOurce processing functions 
source("likelihood.R")
source('get_cox_efficacy.R')
source("format_model_output.R")
source("get_incidence.R")
source("analyse_model_output.R")

# set base parameters 
n_particles = 1L
n_threads = 1L
burnints = 50#90
threshold = 5000
tstep = 1
t_liverstage = 8
country_to_run = orderlyparams$country_to_run
VB = 1e6
divide = if(tstep == 1) 2 else 1

# Set up base inputs (these don't vary across parameter sweep)
base_inputs <- list(
  trial_timesteps = orderlyparams$trial_ts,#365*3,
  burnin = burnints,
  threshold = threshold,
  VB = VB,
  tstep = tstep,
  t_liverstage = t_liverstage,
  country = country_to_run
)

# Create parameter grid
set.seed(123)
# need to make sure to add in a removal or addition of every other day probability of a bite! -- it should be every 2 days 
if (orderlyparams$country_to_run =='BF' | orderlyparams$country_to_run =='Mali'){
  
  A <- randomLHS(n = orderlyparams$n_param_sets, k = 5) # n different sets of parameters, with k parameters to change
  A[,1] <- qunif(A[,1], 2, 20) # min and max values -- this is max kill rate
  A[,2] <- qunif(A[,3], 5, 40) # min and max values -- this is lambda
  A[,3] <- qunif(A[,4], 0.05, 10) # min and max values -- this is kappa
  A[,4] <- round(qunif(A[,4], -50, 50), 0) # lag in days of p of an infectious bite -- this will move the curves to the right 
  A[,5] <- qunif(A[,5], 0.1, 1) # this is the scaling factor for the probability of infectious bite 
  # A
  colnames(A) <- c('max_SMC_kill_rate', 'lambda', 'kappa', 'lag_p_bite', 'p_bite_scaler')
  params_df <- as.data.frame(A)
  params_df$sim_id <- paste0('parameter_set_',rownames(params_df),"_", country_to_run, "_", orderlyparams$sim_allow_superinfections)
  
  # Lag the probability of infectious bites 
  # Get unique lag values from parameter set
  unique_lags <- unique(params_df$lag_p_bite)
  
  # Pre-compute the lagged probabilities for both countries
  # each is a list of the probabilities for each of the n parameter sets in A
  if(country_to_run == 'BF') {
    prob_bite_BFA <- readRDS('prob_bite_BFA.rds') 
    prob_bite_BFA$prob_infectious_bite <- prob_bite_BFA$prob_infectious_bite * params_df$p_bite_scaler
    p_bitevector <-calc_lagged_vectors(prob_bite_BFA, unique_lags, burnints = burnints)
  } else if(country_to_run == 'Mali') {
    prob_bite_MLI <- readRDS('prob_bite_MLI.rds') 
    prob_bite_MLI$prob_infectious_bite <- prob_bite_MLI$prob_infectious_bite * params_df$p_bite_scaler
    p_bitevector <-calc_lagged_vectors(prob_bite_MLI, unique_lags, burnints = burnints) 
  }
  
  params_df$season_start_day <- NA
  
} else {
  param_ranges <- list(
    max_SMC_kill_rate = c(3, 25),      # 
    lambda = c(0.1, 45),   
    kappa = c(0.5, 5),   
    season_start_day = c(20, 120)            # days
  )
  # Generate LHS samples
  A <- randomLHS(orderlyparams$n_param_sets, 4)
  # Scale to parameter ranges
  params_df <- data.frame(
    max_SMC_kill_rate = qunif(A[,1], param_ranges$max_SMC_kill_rate[1], param_ranges$max_SMC_kill_rate[2]),
    lambda = qunif(A[,2], param_ranges$lambda[1], param_ranges$lambda[2]),
    kappa = qunif(A[,3], param_ranges$kappa[1], param_ranges$kappa[1]),
    season_start_day = round(qunif(A[,4], param_ranges$season_start_day[1], param_ranges$season_start_day[2]))
  )
  params_df$sim_id <- paste0('parameter_set_', rownames(params_df),"_", country_to_run, "_", orderlyparams$sim_allow_superinfections)
  
  prob_bite_generic <- readRDS('prob_bite_generic.rds')
  p_bitevector <-calc_lagged_vectors(prob_bite_generic, 0, burnints = burnints) # no lagged values 
  
  params_df$lag_p_bite <- 0
  
  # SMC delivery 
  params_df <- params_df %>%
    rowwise() %>%
    mutate(smc_dose_days = list(c(seq(season_start_day, season_start_day + 120 - 1, 30),
                          seq(season_start_day + 365, season_start_day + 365 + 120 - 1, 30),
                          seq(season_start_day + 365*2, season_start_day  + 365*2 + 120 - 1, 30))))

}


# Add to parameter dataframe by matching lag values (always 0 for generic)
parameters_df <- params_df %>%
  mutate(
    p_bite = map(lag_p_bite, ~p_bitevector[[paste0("lag_", .x)]])
  )

# Make list of parameters instead of df 
params_list <- split(parameters_df, seq(nrow(parameters_df)))

if(orderlyparams$country_to_run != 'generic'){
  # Create a metadata child data frame from trial data 
  metadata_df <- children %>%
    # filter to country we are running 
    filter(country == country_to_run) %>%
    mutate(start_of_fu = as.Date('2017-04-01'),  # this is the start of the cohort simulation. 
           start_to_v3 = v3_date - start_of_fu, # 
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
           smc_dose_days = if_else(arm == 'rtss', list(orderlyparams$trial_ts+burnints), list(smc_dose_days))
    ) %>%
    filter(!is.na(vaccination_day)) %>%
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
    rid = seq(max(metadata_df$rid) + 1, max(metadata_df$rid) + 850),
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
    v1_date = as.Date('2017-04-01') # this is because v1_date is used as start of FU
  ) %>%
    mutate(rid_original = paste0('noint', rid))
  
  metadata_df <- full_join(metadata_df, noint_arm)
} else if (orderlyparams$country_to_run == 'generic'){
  # Day of intervention (0 = start of follow-up; - values are before follow-up; + values after follow-up) - where 0 is also end of burnin
  # these get converted later to the correct directon - i.e. vaccine before follow-up will be +, smc before follow up will be -
  # vax_day is the 3rd primary dose (when we assume that efficacy begins)
  vax_day = -10 # unlike hte model sim, this is in days (not timesteps)
  N = 3200
  
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
    mutate(t_to_boost1 = orderlyparams$trial_ts , 
           t_to_boost2 = orderlyparams$trial_ts + burnints,
           country = country_to_run)

} 
# Run parameter sweep
# results <- run_parameter_sweep(metadata_df, 
#                                parameters_df, 
#                                base_inputs,
#                                output_dir = "simulation_outputs",
#                                parallel = FALSE, # running in parallel doesn't work at hte moment  
#                                allow_superinfections = orderlyparams$sim_allow_superinfections)

# Create folder to save outputs 
output_dir = 'simulation_outputs'
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
# Save parameter grid
saveRDS(parameters_df, file.path(output_dir, "parameter_grid.rds"))
saveRDS(base_inputs, file.path(output_dir, "base_inputs.rds"))
cyphr::encrypt(saveRDS(metadata_df, file.path(output_dir, "metadata_df.rds")), key)

# Set up cluster to run parameter sets in parallel
# Latin Hypercube sampling to do global exploration of parameter space
cluster_cores <- Sys.getenv("CCP_NUMCPUS")
if(cluster_cores == ''){
  cluster_cores = 8
}
message('number of cores: ', cluster_cores)
if (cluster_cores == "") {
  message("running in serial (on a laptop?)")
  message("Running ", nrow(params_df), " simulations sequentially")
  
  results2 <- lapply(params_list, 
                     function(params_row) {
                       run_process(params_row,
                                   metadata_df,
                                   base_inputs,
                                   monthly_inci_trial)
                       # run_cohort_simulation( 
                       #   params_row = params_row,
                       #   metadata_df = metadata_df,
                       #   base_inputs = base_inputs,
                       #   output_dir = output_dir,
                       #   allow_superinfections = orderlyparams$sim_allow_superinfections)
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
    library(orderly2)
    library(retry)
    library(cyphr)
    library(survival)
    library(broom)
    library(survminer)
    
    source('cohort_sim_utils.R')
    source('helper_functions.R')
    source("rtss.R")
    source('run_process.R')
    source("likelihood.R")
    source('get_cox_efficacy.R')
    source("format_model_output.R")
    source("get_incidence.R")
    source("analyse_model_output.R")
    
    TRUE
  })
  
  parallel::clusterExport(cl, c("metadata_df", "base_inputs", "orderlyparams", "gen_bs",
                                "n_particles", "n_threads", "burnints", "threshold", "tstep",
                                "t_liverstage", "country_to_run", "VB", "divide",
                                "monthly_inci_trial"))
  
  results2 <- parallel::clusterApply(cl,
                                       params_list,
                                       function(params_row) {
                                         run_process(params_row,
                                                     metadata_df,
                                                     base_inputs,
                                                     monthly_inci_trial)
                                       }
    )
  # if(orderlyparams$country_to_run == 'smc_fitting'){
  # results2 <- parallel::clusterApply(cl,
  #                                    params_list,
  #                                    function(params_row) {
  #                                      o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
  #                                                                 metadata_df,
  #                                                                 base_inputs,
  #                                                                 output_dir = 'simulation_outputs',
  #                                                                 allow_superinfections = TRUE,
  #                                                                 return_parasitemia = FALSE,
  #                                                                 save_outputs = TRUE)
  # 
  #                                      eff <- calc_efficacy(o$infection_records)
  # 
  #                                      return(list(efficacy = eff,
  #                                                  params = params_row))
  #                                    }
  # )}
  
  parallel::stopCluster(cl)
  
}

# parameters_ll <- unlist(results2)

saveRDS(parameters_ll, file = 'parameters_ll.rds')

parameters_ll <- parameters_ll[!grepl("^p_bite.prob", names(parameters_ll))]
parameters_ll <- parameters_ll[!grepl("^p_bite.date", names(parameters_ll))]
parameters_ll <- parameters_ll[!grepl("^p_bite.year", names(parameters_ll))]
parameters_ll <- parameters_ll[!grepl("^p_bite.rainfall_year", names(parameters_ll))]
library(tidyverse)

# find positions of all sim_id entries
sim_positions <- which(names(parameters_ll) == "sim_id")

# define start and end positions for each parameter group
start_positions <- c(1, sim_positions[-length(sim_positions)] + 1)
end_positions <- sim_positions

# split the flat list into chunks based on these start/end positions
param_groups <- map2(start_positions, end_positions, ~parameters_ll[.x:.y])

# convert each chunk into a one-row tibble and bind together
parameters_df <- map_dfr(param_groups, as_tibble_row)

# convert numeric columns (everything except sim_id) to numeric
parameters_df <- parameters_df %>%
  mutate(across(-sim_id, ~ suppressWarnings(as.numeric(.x))))
