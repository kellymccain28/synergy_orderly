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

orderly_strict_mode()
orderly_parameters(N = NULL, # size of cohort population 
                   trial_ts = NULL)#, # trial timesteps in cohort simulation
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

orderly_dependency(name = 'trial_results',
                   "latest()",
                   c('data/primary_persontime.rds'))

key <- cyphr::data_key()
trial_df <- cyphr::decrypt(readRDS("data/primary_persontime.rds"), key)

# Source antibody function
source("rtss.R")
# Source helper functions
source("helper_functions.R")
# Load the within-host model 
gen_bs <- odin2::odin("smc_rtss.R")
# Source the utils functions 
source("cohort_sim_utils.R")

# Make weights to prevent homogenoeous transmission 
weights <- rexp(N) # Generate random weights to sample children at different probabilities (could also use rpois(N) or rexp(N))
weights <- weights / sum(weights) # normalize to sum to 1 


# Create parameter grid
# max_SMC_kill_rate_vec = c(4,6,8,10)
# SMC_decay_vec = c(0.05, 0.04, 0.03)
# season_start_day_vec = c(20, 50, 80, 110)
# season_length_vec = 30*5

set.seed(123)

A <- randomLHS(n = 3, k = 3) # n different sets of parameters, with k parameters to change
A[,1] <- qunif(A[,1], 4, 10)
A[,2] <- qunif(A[,2], 0.03, 0.1)
A[,3] <- round(qunif(A[,3], 20, 110),0)
A
colnames(A) <- c('max_SMC_kill_rate', 'SMC_decay', 'season_start_day')
params_df <- as.data.frame(A)
# params_df <- expand.grid(
#   max_SMC_kill_rate = max_SMC_kill_rate_vec,
#   SMC_decay = SMC_decay_vec,
#   season_start_day = season_start_day_vec#,
#   # season_length = season_length_vec
# )



# set base parameters 
n_particles = 1L
n_threads = 1L
burnints = 50
threshold = 5000
tstep = 1
t_liverstage = 8
vax_day = -10
country = 'BF'
smc_interval = 30
VB = 1e6
divide = if(tstep == 1) 2 else 1

# Probability of infectious bite 
prob_bite_BFA <- readRDS('prob_bite_BFA.rds')
prob_bite_MLI <- readRDS('prob_bite_MLI.rds')
# add in median value of probability of infectious bite at the beginning of the vector for the burnin period
p_bite_bfa <- c(rep(median(prob_bite_BFA$prob_infectious_bite), burnints), prob_bite_BFA$prob_infectious_bite)
p_bite_mli <- c(rep(median(prob_bite_MLI$prob_infectious_bite), burnints), prob_bite_MLI$prob_infectious_bite)


# Set up base inputs (these don't vary across parameter sweep)
base_inputs <- list(
  N_children = N,
  trial_timesteps = trial_ts,#365*3,
  burnin = burnints,
  p_bite_bfa = p_bite_bfa,  # your loaded data
  p_bite_mli = p_bite_mli,  # your loaded data
  threshold = threshold,
  VB = VB,
  tstep = tstep,
  t_liverstage = t_liverstage,
  weights = weights,        # your generated weights
  vax_day = vax_day,
  country = country,
  smc_interval = smc_interval
)


# Run parameter sweep
results <- run_parameter_sweep(params_df, 
                               base_inputs,
                               output_dir = "simulation_outputs",
                               parallel = FALSE # running in parallel doesn't work at hte moment  
                               )

