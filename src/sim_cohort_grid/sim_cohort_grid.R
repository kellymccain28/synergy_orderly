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
orderlyparams <- orderly_parameters(N = NULL, # size of cohort population 
                   trial_ts = NULL,
                   sim_allow_superinfections = NULL)#, # trial timesteps in cohort simulation
                   # burnin = NULL, # number of days in cohort simulation to use for burnin
                   # max_SMC_kill_rate = NULL,
                   # SMC_decay = NULL,
                   # season_start_day = NULL,
                   # season_length = NULL,
                   # smc_interval = NULL)

orderly_shared_resource("smc.R",
                        "rtss.R",
                        "smc_rtss.R",
                        "helper_functions.R")

orderly_resource("cohort_sim_utils.R")

orderly_dependency(name = 'fit_rainfall',
                   "latest()",
                   c("prob_bite_BFA.rds",
                     "prob_bite_MLI.rds"))

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

# Make weights to prevent homogenoeous transmission 
weights <- rexp(orderlyparams$N) # Generate random weights to sample children at different probabilities (could also use rpois(N) or rexp(N))
weights <- weights / sum(weights) # normalize to sum to 1 

# set base parameters 
n_particles = 1L
n_threads = 1L
burnints = 50
threshold = 5000
tstep = 1
t_liverstage = 8
vax_day = -10
country_to_run = 'BF'
smc_interval = 30
VB = 1e6
divide = if(tstep == 1) 2 else 1


# Set up base inputs (these don't vary across parameter sweep)
base_inputs <- list(
  N_children = orderlyparams$N,
  trial_timesteps = orderlyparams$trial_ts,#365*3,
  burnin = burnints,
  threshold = threshold,
  VB = VB,
  tstep = tstep,
  t_liverstage = t_liverstage,
  weights = weights,        
  vax_day = vax_day,
  country = country_to_run,
  smc_interval = smc_interval
)

# Create parameter grid
set.seed(123)

A <- randomLHS(n = 20, k = 4) # n different sets of parameters, with k parameters to change
A[,1] <- qunif(A[,1], 2, 15) # min and max values -- this is max kill rate
A[,2] <- qunif(A[,2], 0.01, 0.1) # min and max values -- this is SMC decay rate
A[,3] <- 122#round(qunif(A[,3], 90, 150),0) # min and max values -- this is day of season start, which is also the first day of SMC -- relative to the first day of FU which is 1 April 2017
A[,4] <- round(qunif(A[,4],0, 60), 0) # lag in days of p of an infectious bite -- this will move the curves to the right 
# A
colnames(A) <- c('max_SMC_kill_rate', 'SMC_decay', 'season_start_day', 'lag_p_bite')
params_df <- as.data.frame(A)

# Lag the probability of infectious bites 
prob_bite_BFA <- readRDS('prob_bite_BFA.rds')
prob_bite_MLI <- readRDS('prob_bite_MLI.rds')
# Get unique lag values from parameter set
unique_lags <- unique(params_df$lag_p_bite)

#Calculate lagged vectors for all unique lags
calc_lagged_vectors <- function(prob_data, lags, start_date = '2017-04-01', 
                                      end_date = '2020-04-01', burnints) {
  
  lag_list <- map(lags, function(lag_val) {
    prob_lagged <- prob_data %>% 
      mutate(prob_lagged = lag(prob_infectious_bite, n = lag_val),
             date_lagged = lag(date, n = lag_val))
    
    prob_filtered <- prob_lagged[prob_lagged$date_lagged >= start_date & 
                                   prob_lagged$date_lagged < end_date & 
                                   !is.na(prob_lagged$date_lagged),]
    
    c(rep(median(prob_filtered$prob_lagged, na.rm = TRUE), burnints), 
      prob_filtered$prob_lagged)
  })
  
  names(lag_list) <- paste0("lag_", lags)
  return(lag_list)
}

# Pre-compute the lagged probabilities for both countries
bfa_vectors <- calc_lagged_vectors(prob_bite_BFA, unique_lags, burnints = burnints)
mli_vectors <- calc_lagged_vectors(prob_bite_MLI, unique_lags, burnints = burnints)

# Add to parameter dataframe by matching lag values
parameters_df <- params_df %>%
  mutate(
    p_bite_bfa = map(lag_p_bite, ~bfa_vectors[[paste0("lag_", .x)]]),
    p_bite_mli = map(lag_p_bite, ~mli_vectors[[paste0("lag_", .x)]])
  )

# Create a metadata child data frame from trial data 

metadata_df <- children %>%
  mutate(start_of_fu = v1_date, 
         start_to_v3 = v3_date - start_of_fu, # then the negative version of this is what we need to go into the vax_day var (or could change it to be pos always)
         # vaccination day is the day of vaccination relative to the start of follow-up 
         vaccination_day = as.numeric(-start_to_v3),
         PEV = ifelse(arm == 'smc', 0, 1),
         SMC = ifelse(arm == 'rtss', 0, 1),
         t_to_boost1 = as.numeric(boost1_date - v3_date),
         t_to_boost2 = as.numeric(boost2_date - v3_date)) %>%
  # Calculate timings of smc for each year
  mutate() %>%
  select(-v3_date, -v1_date) %>%
  filter(!is.na(vaccination_day)) %>%
  # for second booster, say if it is missing then second booster is much later so that it is after the follow-up time is over
  mutate(t_to_boost2 = ifelse(is.na(t_to_boost2), 1400, t_to_boost2),
         t_to_boost2 = ifelse(is.na(t_to_boost2), 1500, t_to_boost2)) %>%
  # update RIDs to be numeric 
  mutate(rid_original = rid) %>%
  group_by(country) %>%
  mutate(rid = row_number())



# Run parameter sweep
results <- run_parameter_sweep(metadata_df, 
                               parameters_df, 
                               base_inputs,
                               output_dir = "simulation_outputs",
                               parallel = FALSE, # running in parallel doesn't work at hte moment  
                               allow_superinfections = orderlyparams$sim_allow_superinfections)

