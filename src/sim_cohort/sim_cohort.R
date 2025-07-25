# Script to simulate a cohort 
# 17 Feb 2025
# Kelly McCain
# Last updated: 7 July 2025

library(odin2)
library(ggplot2)
library(dust2)
library(tidyverse)
library(mgcv)
library(umbrella)

orderly_strict_mode()
orderly_parameters(N = NULL, # size of cohort population 
                   trial_ts = NULL, # trial timesteps in cohort simulation 
                   burnin = NULL, # number of days in cohort simulation to use for burnin
                   max_SMC_kill_rate = NULL,
                   SMC_decay = NULL,
                   season_start_day = NULL,
                   season_length = NULL,
                   smc_interval = NULL)

# orderly_artefact(description = 'input and output datasets and diagnostic plots for cohort simulation',
#                  files = c(
#                    'inputs.rds',
#                    'outputs/metadata_children.rds',
#                    'outputs/parasitemia.rds',
#                    'outputs/infection_records.rds',
#                    "outputs/plots/parasitemia_random.png",
#                    'outputs/plots/smc_timings.png',
#                    'outputs/plots/proportion_detectable.png',
#                    "outputs/plots/incidence.png",
#                    'outputs/plots/cum_infections.png'
#                    ))


orderly_shared_resource("smc.R",
                        "rtss.R",
                        "smc_rtss.R",
                        "helper_functions.R")

orderly_dependency(name = 'fit_rainfall',
                   "latest()",
                   c("prob_bite_BFA.rds",
                     "prob_bite_MLI.rds"))

# Source antibody function
source("rtss.R")
# Source helper functions
source("helper_functions.R")
# Load the within-host model 
gen_bs <- odin2::odin("smc_rtss.R")

# Create folders to save outputs 
if (!dir.exists('outputs')) {dir.create("outputs")}
if (!dir.exists('outputs/plots')) {dir.create("outputs/plots")}


# Cohort information 
# N <- 5000
trial_ts <- 365*3 # Total timesteps in the cohort simulation
burnin = 50 # number of simulation timesteps before the interventions are delivered and we start follow-up

# fourier series for seasonality 
# rainfall_coefs <- c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919) # seasonal
# rainfall <- umbrella:::fourier_predict(coef = rainfall_coefs,
#                             t = seq(1, 365),
#                             floor = 0.001)
# rainfall$profile <- rainfall$profile / 150
# p_bite <- rep(rainfall$profile, (trial_ts)/365 + 1)# 0.002 # probability of infectious bite over 1 day 
# plot(p_bite[1:(365*3)])

# Probability of infectious bite 
prob_bite_BFA <- readRDS('prob_bite_BFA.rds')
prob_bite_MLI <- readRDS('prob_bite_MLI.rds')
# add in median value of probability of infectious bite at the beginning of the vector for the burnin period
p_bite_bfa <- c(rep(median(prob_bite_BFA$prob_infectious_bite), burnin), prob_bite_BFA$prob_infectious_bite)
p_bite_mli <- c(rep(median(prob_bite_MLI$prob_infectious_bite), burnin), prob_bite_MLI$prob_infectious_bite)

# abline(v = 50, col = 'red')
threshold <- 5000#100 # PB threshold per microL for when this is a detectable infection 

n_particles = 1L
n_threads = 1L
# t_inf = 2 * 365
# ts = 300 # number of days in infection model
VB = 1e6 # blood volume in microL
det_mode = FALSE # deterministic mode
tstep <- 1 # timestep 1 or 2 days -- ie. is the model itself ocunting every 2 days or is it by 1 day and then we mutliply by 2 after?
divide <- if(tstep == 1) 2 else 1
tt <- seq(0, ts/divide, by = tstep)# timesteps that the individual within-host model will run (depending on timestep)
# max_SMC_kill_rate = 10#6
# SMC_decay = 0.05
t_liverstage <- 8 #/ (2 / tstep) # time in days spent in liver stage / length of each timetsep in days 

# SMC parameters
# season_start_day <- 213 # August 1 
# season_length <- 120
# smc_interval <- 30

weights <- rexp(N) # Generate random weights to sample children at different probabilities (could also use rpois(N) or rexp(N))
weights <- weights / sum(weights) # normalize to sum to 1 
# ggplot() + geom_density(aes(x = weights))

# Make a matrix of potential values for each parameter 
# max_SMC_kill_rate = c(4,6,8,10)
# SMC_decay = c(0.05, 0.03)
# season_start_day = c(50, 80, 110)
# season_length = 30*5
# smc_interval = 30

# params <- expand.grid(
#   max_SMC_kill_rate = max_SMC_kill_rate,
#   SMC_decay = SMC_decay,
#   season_start_day = season_start_day,
#   season_length = season_length, 
#   smc_interval
# )

# Make dataframe of the different inputs to the model 
inputs <- data.frame(
  N_children = N, 
  trial_timesteps = I(list(trial_ts)), # this I() function prevents automatic unlisting 
  rainfall = I(list(rainfall)),
  p_bite = I(list(p_bite)),
  threshold = threshold, 
  VB = VB, 
  tstep = tstep,
  t_liverstage = t_liverstage,
  burnin = burnin, 
  max_SMC_kill_rate = if (length(max_SMC_kill_rate) == 1) max_SMC_kill_rate else I(list(max_SMC_kill_rate)),
  SMC_decay = if (length(SMC_decay) == 1) SMC_decay else I(list(SMC_decay)), 
  season_start_day = season_start_day,
  season_length = season_length, 
  smc_interval = smc_interval, 
  weights = I(list(weights))
)

saveRDS(inputs, file = 'inputs.rds')


# Run trial simulation where we append infections at the end of the df ----
# Initialize storage
infection_records <- data.frame(
  child_id = integer(),             # individual id of each child
  BSinfection_day = integer(),      # day of the blood stage infection (simulation day)
  threshold_day = integer(),        # day since blood-stage begins that infection is detected
  t_toreach_threshold = integer(),  # day since infectious bite that infection is detected
  vaccination_day = integer(),       # Store vaccination day
  PEV = integer(), 
  SMC = integer(), 
  intervention = character(),
  p_bite = double()
)

# Day of intervention (0 = start of follow-up; - values are before follow-up; + values after follow-up) - where 0 is also end of burnin
# these get converted later to the correct directon - i.e. vaccine before follow-up will be +, smc before follow up will be -
# vax_day is the 3rd primary dose (when we assume that efficacy begins)
vax_day = -10 # unlike hte model sim, this is in days (not timesteps)

metadata_child <- data.frame(
  child_id = 1:N,
  vaccination_day = vax_day,#sample(-1:0, N, replace = TRUE), the model takes abs(days before start of follow-up)
  PEV = rbinom(N, 1, 0.5),#rep(0,N),#
  SMC = c(rep(1, N/2), rep(0, N/2))#rbinom(N, 1, 0.5)
) %>%
  mutate(intervention = case_when(
    PEV == 1 & SMC == 1 ~ 'vaxsmc',
    PEV == 1 & SMC == 0 ~ 'vax',
    PEV == 0 & SMC == 1 ~ 'smc', 
    TRUE ~ 'none'))

saveRDS(metadata_child, 'outputs/metadata_children.rds')#'smctime_',smc_day, 


# Initialize list to store parasitaemia trajectories per child
parasitemia_storage <- vector("list", length = trial_ts + burnin)
country = 'BF'
for (t in 1:(trial_ts + burnin)){
  
  p_bite <- if(country == 'BF') p_bite_bfa else p_bite_mli
  
  # Determine number of new infectious bites on time t
  n_infectious_bites <- rpois(1, lambda = N * p_bite[t]) # or could do rbinom(1, N, p_bite)
  
  if (n_infectious_bites == 0) next
  
  # Sample random children to be bitten 
  bites <- sample(1:N, size = n_infectious_bites, prob = weights, replace = TRUE) #prob = weights,
  message('There are ', length(bites), " mosquito bites on time ", t, ' out of ', (trial_ts+burnin))
  
  bit_kids <- unique(bites)
  message('There are ', length(bit_kids), " unique kids bitten on time ", t, ' out of ', (trial_ts+burnin))
  
  if(length(bit_kids) > 0) {
    # Estimate antibody level for each bitten child 
    
    if(t < burnin){
      # No interventions during burnin
      PEV_vec <- rep(0, length(bit_kids))
      SMC_vec <- rep(0, length(bit_kids))
      # these both need a placeholder value (won't matter if PEV = 0 or SMC = 0)
      t_since_vax_vec <- rep(0, length(bit_kids))
      
      infection_start_day = t
    } else {# Find time in days between vaccination and infection
      
      # Vectorized intervention calculations
      kid_metadata <- metadata_child[metadata_child$child_id %in% bit_kids, ]
      
      # Create named vectors
      pev_lookup <- setNames(kid_metadata$PEV, kid_metadata$child_id)
      smc_lookup <- setNames(kid_metadata$SMC, kid_metadata$child_id)
      vax_day_lookup <- setNames(kid_metadata$vaccination_day, kid_metadata$child_id)
      
      # Map to bit_kids order
      PEV_vec <- pev_lookup[as.character(bit_kids)]
      SMC_vec <- smc_lookup[as.character(bit_kids)]
      t_since_vax_vec <- (t - burnin) - vax_day_lookup[as.character(bit_kids)]
      
      # Get intervention status
      message('PEV:', PEV_vec)
      message('SMC:', SMC_vec)
    }
    # message('pev=', PEV, " and smc=", SMC)
    # Run the within-host simulation for the kids who have been bitten
    params_df <- data.frame(
      PEV_on = PEV_vec,
      SMC_on = SMC_vec,
      t_inf = t_since_vax_vec,
      # smc_timing = smc_time_vec,
      infection_start_day = t,
      season_start_day = season_start_day,
      season_length = season_length, 
      smc_interval = smc_interval,
      child_id = bit_kids
    )
    # Find time length until end of follow-up to run infection sim  
    tt_until_end_cohort <- seq(0, ((trial_ts + burnin) - t)/divide, by = tstep)#
    # message('time until end of cohort, run infection sim for: ', tt_until_end_cohort, ' timesteps')
 
    outputs <- pmap(params_df, 
                    function(PEV_on, 
                             SMC_on, 
                             t_inf, 
                             infection_start_day,
                             season_start_day, 
                             season_length, 
                             smc_interval,
                             child_id) {
                      result <- run_process_model(
        PEV_on = PEV_on,
        SMC_on = SMC_on,
        t_inf = t_inf,
        VB = VB,
        tt = tt_until_end_cohort,
        max_SMC_kill_rate = max_SMC_kill_rate,
        SMC_decay = SMC_decay,
        # smc_timing = smc_timing,
        infection_start_day = infection_start_day, # external time that infection begins 
        season_start_day = season_start_day, # day of season start relative to Jan 1 of external year 
        season_length = season_length, # ~4 months (120 days)
        smc_interval= smc_interval # how often are SMC rounds (days)
      )
                      
      result$trajectory$child_id <- child_id
      
      return(result)
    })
    
    # Vectorized data frame creation for infection records
    new_records <- data.frame(
      time_ext = t - burnin,                                                      # external cohort time (t is the cohort time, then we wnat to scale to be + if after burnin)
      child_id = bit_kids,                                                        # number of bit kids
      infectious_bite_day = (t - burnin) - t_liverstage,                          # bitten on day t, then assuming the liver stage takes t_liverstage days
      BSinfection_day = (t - burnin),                                             # after liver stage, the BS begins #+ t_liverstage
      threshold_day = sapply(outputs, function(x) x$threshold_day),               # days since BS starts that threshold is reached
      detection_day = (t - burnin) + sapply(outputs, function(x) x$threshold_day),# day in cohort simulation that threshold is reached, threshold is the day since BS infection that reaches threshold #+ t_liverstage 
      t_toreach_threshold = sapply(outputs, function(x) x$threshold_day) + t_liverstage,   # time to reach threshold value / detection since the bite  
      vaccination_day = if(t < burnin) rep(vax_day, length(bit_kids)) else kid_metadata$vaccination_day,   # day of vaccination relative to the start of follow-up (day 0 external time)
      prob_bite = p_bite[t]
    ) 
    
    infection_records <- rbind(infection_records, new_records)
    
    # Vectorized parasitemia storage creation
    parasitemia_data <- map2(outputs, bit_kids, function(output, kid) {
      
      output$trajectory %>%
        mutate(
          day1_BSinfection = t - burnin,#+ t_liverstage 
          detection_day = t + (output$threshold_day) - burnin,# do not need to multiply threshold day by 2 here since already done in run_process_model #+ t_liverstage 
          time_ext = time*2 + (t - burnin),# + infection_start_day/2,# external time should be dependent on when infection was, relative to external time(inf_start_day aka t); time (model time) has already been multiplied by 2 and related to infection time  #time*2 + (t - 1) - burnin,#+ t_liverstage
          intervention = metadata_child[metadata_child$child_id == kid, ]$intervention,
          t = t
        )
    })
    
    # Store all parasitemia data for this time step
    parasitemia_storage[[t]] <- parasitemia_data
    
  }
}

# Flatten parasitemia storage at the end (more efficient than appending)
parasitemia_storage <- unlist(parasitemia_storage, recursive = FALSE)

parasitemia_df <- bind_rows(parasitemia_storage) %>%
  filter(detection_day > 0 | is.na(detection_day)) %>%
  group_by(child_id, day1_BSinfection) %>%
  mutate(child_dayinf = paste0(child_id, ", day ", day1_BSinfection),
  ) %>%
  mutate(det = ifelse(!is.na(detection_day), 1, 0)) #%>% # this is the day of the cohort simulation that the infection is found --

saveRDS(parasitemia_df, file = 'outputs/parasitemia.rds')


# Add in relevant metadata to infection records df
child_counts <- metadata_child %>%
  distinct(child_id, intervention) %>%
  count(intervention, name = "children_in_group")

infection_records2 <- infection_records %>% 
  filter(detection_day > 0 | is.na(detection_day)) %>% # this is the day of the cohort simulation that the infection is found --
  left_join(metadata_child %>% select(-c(vaccination_day)), by = c('child_id')) %>%
  left_join(child_counts, by = 'intervention') %>%
  mutate(detectable = ifelse(is.na(threshold_day), 0, 1)) %>%
  group_by(intervention) %>%
  arrange(detection_day) %>%
  mutate(cumul_inf = cumsum(detectable))

saveRDS(infection_records2, file = 'outputs/infection_records.rds')


# Make plots 

# Parasitemia over time by child in sample of 100 random children per group
parasitemia_random <- ggplot(parasitemia_df[parasitemia_df$child_id %in% (parasitemia_df %>% distinct(child_id, intervention) %>% group_by(intervention) %>% slice_sample(n=50) %>% pull(child_id)), ]) +
  geom_line(aes(x=time_ext/365.25, y = parasites, group = as.factor(child_dayinf ), color = as.factor(det )), alpha = 0.5) +
  # geom_line(aes(x=time_ext/365.25, y = exp(SMC_kill_rateout), group = as.factor(child_dayinf ), color = as.factor(det )), alpha = 0.5) +
  scale_y_log10() +
  geom_hline(aes(yintercept = 10), linetype = 2, color = 'darkred', linewidth = 1) + # this is the detection limit (following Challenger et al.)
  geom_hline(aes(yintercept = 1e-5), linetype = 2, color = 'darkgreen', linewidth = 1) + # this is the clearance threshold
  geom_vline(aes(xintercept = 0), color = 'darkturquoise', linetype = 2, linewidth = 1) +
  facet_wrap(~intervention) + 
  theme_bw() +
  theme(legend.position = 'none')

#check SMC timings  
smc_timings <- ggplot(parasitemia_df[parasitemia_df$child_dayinf %in% (parasitemia_df %>% group_by(intervention) %>% slice_sample(n=3) %>% pull(child_dayinf)) &
                        parasitemia_df$intervention %in% c('smc','vaxsmc'),]) +#sample(parasitemia_df$child, size = 50)
  # geom_line(aes(x = time_ext/365.25, y = log10(numkillSMC), color = as.factor(detection_day))) +
  geom_line(aes(x = time_ext/365.25, y = prob_smckill, color = as.factor(detection_day))) +
  geom_line(aes(x = time_ext/365.25, y = log10(parasites)/2), color = 'black') +
  geom_line(aes(x = time_ext/365.25, y = growth/2), color = 'pink', alpha = 0.5) +
  geom_line(data = as.data.frame(rainfall), aes(x = t/365.25, y = profile*100)) +
  theme_bw() +
  facet_wrap(~child_dayinf + intervention) 


prop_det <- ggplot(infection_records2 ) + 
  geom_bar(aes(x = intervention, group = as.factor(detectable), fill = as.factor(detectable)),
           position = 'fill') + 
  scale_fill_manual(values = c('darkmagenta','goldenrod')) +
  labs(#title = 'Detectable infections in each intervention group',
    fill = 'Detectable')

# incidence
incidence <- ggplot(infection_records2 )+#%>% filter(detection_day <= ts))+
  geom_line(aes(x = detection_day, color=intervention), stat = 'count') + 
  facet_wrap(~intervention) + 
  labs(y = 'N infections',
       x = 'Day since start of follow up period',
       caption = 'Assuming that the liver stage lasts 8 days') + 
  theme(legend.position = 'none')

# Cumulative number of infections 
cum_infections <- ggplot(infection_records2 )+#%>% filter(detection_day <= ts))+
  geom_line(aes(x = detection_day, y = cumul_inf, color=intervention)) + 
  labs(y = 'N infections',
       x = 'Day since start of follow up period',
       caption = 'Assuming that the liver stage lasts 8 days') 

ggsave(filename = "outputs/plots/parasitemia_random.png", parasitemia_random, height = 8, width = 12)
ggsave(filename = 'outputs/plots/smc_timings.png', smc_timings, height = 8, width = 12)
ggsave('outputs/plots/proportion_detectable.png', prop_det, height = 8, width = 12)
ggsave("outputs/plots/incidence.png", incidence, height = 8, width = 12)
ggsave('outputs/plots/cum_infections.png', cum_infections, height = 8, width = 12)

# dev.off()
