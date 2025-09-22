library(odin2)
library(ggplot2)
library(dust2)
library(tidyverse)
library(mgcv)
library(umbrella)
library(future)
library(future.apply)
library(lhs)
library(orderly2)
library(cyphr)

key <- cyphr::data_key()

orderly_strict_mode()
orderlyparams <- orderly_parameters(#N = NULL, # size of cohort population 
                   trial_ts = NULL,# trial timesteps in cohort simulation (inte)
                   sim_allow_superinfections = NULL, # TRUE or FALSE 
                   country_to_run = NULL, # BF or Mali, or if generic, then 'generic'
                   n_param_sets = NULL)

                   

orderly_shared_resource("smc.R",
                        "rtss.R",
                        "smc_rtss.R",
                        "helper_functions.R")

orderly_resource("cohort_sim_utils.R")

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

# Source antibody function
source("rtss.R")
# Source helper functions
source("helper_functions.R")
# Load the within-host model 
gen_bs <- odin2::odin("smc_rtss.R")
# Source the utils functions 
source("cohort_sim_utils.R")

# set base parameters 
n_particles = 1L
n_threads = 1L
burnints = 50
threshold = 5000
tstep = 1
t_liverstage = 8
country_to_run = orderlyparams$country_to_run
VB = 1e6
divide = if(tstep == 1) 2 else 1

# Set up base inputs (these don't vary across parameter sweep)
base_inputs <- list(
  # N_children = nrow(children),
  trial_timesteps = orderlyparams$trial_ts,#365*3,
  burnin = burnints,
  threshold = threshold,
  VB = VB,
  tstep = tstep,
  t_liverstage = t_liverstage,
  weights = weights,        
  country = country_to_run
)

# Create parameter grid
set.seed(123)

A <- randomLHS(n = orderlyparams$n_param_sets, k = 4) # n different sets of parameters, with k parameters to change
A[,1] <- qunif(A[,1], 2, 15) # min and max values -- this is max kill rate
A[,2] <- qunif(A[,3], 10, 30) # min and max values -- this is lambda
A[,3] <- qunif(A[,4], 0.05, 0.3) # min and max values -- this is kappa
A[,4] <- round(qunif(A[,4],0, 100), 0) # lag in days of p of an infectious bite -- this will move the curves to the right 
# A
colnames(A) <- c('max_SMC_kill_rate', 'lambda', 'kappa', 'lag_p_bite')
params_df <- as.data.frame(A)
params_df$sim_id <- paste0('sim_', rownames(params_df))

# Lag the probability of infectious bites 
prob_bite_BFA <- readRDS('prob_bite_BFA.rds')
prob_bite_MLI <- readRDS('prob_bite_MLI.rds')
prob_bite_generic <- readRDS('prob_bite_generic.rds')

# Get unique lag values from parameter set
unique_lags <- unique(params_df$lag_p_bite)

#Calculate lagged vectors for all unique lags
calc_lagged_vectors <- function(prob_data, lags, start_date = as.Date('2017-04-01'), 
                                      end_date = '2020-04-01', burnints) {
  
  lag_list <- map(lags, function(lag_val) {
    prob_lagged <- prob_data %>% 
      mutate(prob_lagged = dplyr::lag(prob_infectious_bite, n = lag_val),
             date_lagged = dplyr::lag(date, n = lag_val))
    
    # Get start date minus burnin 
    start_date_pbite <- start_date - burnints
    prob_filtered <- prob_lagged[prob_lagged$date_lagged >= start_date_pbite & 
                                   prob_lagged$date_lagged < end_date & 
                                   !is.na(prob_lagged$date_lagged),]
    
    # isntead of median, am now filtering to start date - burnin above
    # c(rep(median(prob_filtered$prob_lagged, na.rm = TRUE), burnints), # this is to have a probability of bite before the burnin 
    #   prob_filtered$prob_lagged)
  })
  
  names(lag_list) <- paste0("lag_", lags)
  return(lag_list)
}

# Pre-compute the lagged probabilities for both countries
# each is a list of the probabilities for each of the n parameter sets in A
bfa_vectors <- calc_lagged_vectors(prob_bite_BFA, unique_lags, burnints = burnints)
mli_vectors <- calc_lagged_vectors(prob_bite_MLI, unique_lags, burnints = burnints)
generic_vectors <- calc_lagged_vectors(prob_bite_generic, 0, burnints = burnints)

p_bitevector <- if(country_to_run == 'BF') {
  bfa_vectors
} else if(country_to_run == 'Mali') {
  mli_vectors 
} else if(country_to_run == 'generic') {
  generic_vectors
}

# Add to parameter dataframe by matching lag values
parameters_df <- params_df %>%
  mutate(
    p_bite = map(lag_p_bite, ~p_bitevector[[paste0("lag_", .x)]])#,
    # p_bite_mli = map(lag_p_bite, ~mli_vectors[[paste0("lag_", .x)]])
  )

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
  mutate(rid = as.numeric(str_extract(rid, "[0-9]+"))) 


# Run parameter sweep
results <- run_parameter_sweep(metadata_df, 
                               parameters_df, 
                               base_inputs,
                               output_dir = "simulation_outputs",
                               parallel = FALSE, # running in parallel doesn't work at hte moment  
                               allow_superinfections = orderlyparams$sim_allow_superinfections)

