# Script to simulate a cohort 
# 17 Feb 2025
# Kelly McCain
# Last updated: 26 Sep 2025

library(odin2)
library(ggplot2)
library(dust2)
library(tidyverse)
library(mgcv)
library(umbrella)
library(orderly2)

orderly_strict_mode()
cohortparams <- orderly_parameters(N = NULL, # size of cohort population 
                   trial_ts = NULL, # trial timesteps in cohort simulation 
                   burnin = NULL, # number of days in cohort simulation to use for burnin
                   max_SMC_kill_rate = NULL,
                   smc_lambda = NULL,  # default could be 17
                   smc_kappa = NULL, # default could be 0.28
                   season_start_day = NULL,
                   sim_allow_superinfections = NULL
                   # season_length = NULL,
                   # smc_interval = NULL
                   )


orderly_shared_resource("smc.R",
                        "rtss.R",
                        "smc_rtss.R",
                        "helper_functions.R")

orderly_dependency(name = 'fit_rainfall',
                   "latest()",
                   c("prob_bite_BFA.rds",
                     "prob_bite_MLI.rds",
                     "prob_bite_generic.rds"))

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
# trial_ts <- 365*3 # Total timesteps in the cohort simulation
# burnin = 50 # number of simulation timesteps before the interventions are delivered and we start follow-up

# Probability of infectious bite 
prob_bite_BFA <- readRDS('prob_bite_BFA.rds')
prob_bite_MLI <- readRDS('prob_bite_MLI.rds')
prob_bite_generic <- readRDS('prob_bite_generic.rds')
# add in median value of probability of infectious bite at the beginning of the vector for the burnin period
# p_bite_bfa <- c(rep(median(prob_bite_BFA$prob_infectious_bite), burnin),
#                 prob_bite_BFA$prob_infectious_bite)
# p_bite_mli <- c(rep(median(prob_bite_MLI$prob_infectious_bite), burnin), 
#                 prob_bite_MLI$prob_infectious_bite)
start_cohort <- as.Date('2017-04-01')
start_date_pbite <- start_cohort - cohortparams$burnin
end_date <- start_date_pbite + cohortparams$trial_ts + cohortparams$burnin
p_bite_generic <- prob_bite_generic[prob_bite_generic$date >= start_date_pbite & 
                                      prob_bite_generic$date < end_date & 
                                                 !is.na(prob_bite_generic$date),]

    

threshold <- 5000#100 # PB threshold per microL for when this is a detectable infection 

n_particles = 1L
n_threads = 1L
# t_inf = 2 * 365
ts = 300 # number of days in infection model
VB = 1e6 # blood volume in microL
det_mode = FALSE # deterministic mode
tstep <- 1 # timestep 1 or 2 days -- ie. is the model itself ocunting every 2 days or is it by 1 day and then we mutliply by 2 after?
divide <- if(tstep == 1) 2 else 1
tt <- seq(1, ts/divide, by = tstep)# timesteps that the individual within-host model will run (depending on timestep)
# max_SMC_kill_rate = 10#6
# SMC_decay = 0.05
t_liverstage <- 8 #/ (2 / tstep) # time in days spent in liver stage / length of each timetsep in days 
allow_superinfections = cohortparams$sim_allow_superinfections

weights <- rexp(cohortparams$N) # Generate random weights to sample children at different probabilities (could also use rpois(N) or rexp(N))
weights <- weights / sum(weights) # normalize to sum to 1 
# ggplot() + geom_density(aes(x = weights))


# Make function to do cohort simulations over the grid of parameter values 

# Run trial simulation where we append infections at the end of the df ----
# Initialize storage
infection_records <- data.frame(
  rid = integer(),             # individual id of each child
  BSinfection_day = integer(),      # day of the blood stage infection (simulation day)
  threshold_day = integer(),        # day since blood-stage begins that infection is detected
  t_toreach_threshold = integer(),  # day since infectious bite that infection is detected
  vaccination_day = integer(),       # Store vaccination day
  PEV = integer(), 
  SMC = integer(), 
  arm = character(),
  p_bite = double(),
  recovery_day = double()
)

# Day of intervention (0 = start of follow-up; - values are before follow-up; + values after follow-up) - where 0 is also end of burnin
# these get converted later to the correct directon - i.e. vaccine before follow-up will be +, smc before follow up will be -
# vax_day is the 3rd primary dose (when we assume that efficacy begins)
vax_day = -10 # unlike hte model sim, this is in days (not timesteps)

metadata_child <- data.frame(
  rid = 1:cohortparams$N,
  vaccination_day = vax_day,#sample(-1:0, N, replace = TRUE), the model takes abs(days before start of follow-up)
  PEV = rbinom(cohortparams$N, 1, 0.5),#rep(0,N),#
  SMC = c(rep(1, cohortparams$N/2), rep(0, cohortparams$N/2))#rbinom(N, 1, 0.5)
) %>%
  mutate(intervention = case_when(
    PEV == 1 & SMC == 1 ~ 'vaxsmc',
    PEV == 1 & SMC == 0 ~ 'vax',
    PEV == 0 & SMC == 1 ~ 'smc', 
    TRUE ~ 'none'))

saveRDS(metadata_child, 'outputs/metadata_children.rds')#'smctime_',smc_day, 

# add in generic smc dose days that are evenly spaced throughout the season -- 
# make sure season is aligned with the pbite vector 
season_start_day = cohortparams$season_start_day
smc_delivery_dates <- c(seq(season_start_day, season_start_day + 120 - 1, 30),
                        seq(season_start_day + 365, season_start_day + 365 + 120 - 1, 30),
                        seq(season_start_day + 365*2, season_start_day  + 365*2 + 120 - 1, 30))
dates <- seq(as.Date('2017-02-10'), as.Date('2020-03-31'))
day_since_start <- time_length(interval(as.Date('2017-04-01'), dates), unit = 'days')

calculate_days_since_smc <- function(day_since_start, smc_delivery_dates) {
  # Find the most recent SMC delivery date that is <= current day
  recent_smc <- max(smc_delivery_dates[smc_delivery_dates <= day_since_start], na.rm = TRUE)
  # If no SMC delivery has occurred yet, return NA
  if (is.infinite(recent_smc)) {return(NA)}
  # Return days since most recent SMC delivery
  return(day_since_start - recent_smc)
}

time_since_smc <- sapply(day_since_start, calculate_days_since_smc, smc_delivery_dates)
smckillvec <- cohortparams$max_SMC_kill_rate * exp(-(time_since_smc/ cohortparams$smc_lambda)^cohortparams$smc_kappa)  # calculate kill rate with hill function
smckillvec <- ifelse(is.na(smckillvec), 0, smckillvec)                                       # change NAs (when there is no SMC) to 0 
smckillvec <- c(rep(0, cohortparams$burnin), smckillvec) 
smckillvec <- smckillvec[seq_along(smckillvec) %% 2 == 1]  
metadata_child$smckillvec <- I(list(smckillvec))

# Make dataframe of the different inputs to the model 
inputs <- list(
  N_children = cohortparams$N, 
  trial_timesteps = cohortparams$trial_ts,
  p_bite = I(list(p_bite_generic$prob_infectious_bite)),
  threshold = threshold, 
  VB = VB, 
  tstep = tstep,
  t_liverstage = t_liverstage,
  burnin = cohortparams$burnin, 
  max_SMC_kill_rate = I(list(cohortparams$max_SMC_kill_rate)),
  smc_lambda = cohortparams$smc_lambda,
  kappa = cohortparams$smc_kappa,
  season_start_day = cohortparams$season_start_day,
  weights = I(list(weights)),
  smckillvec = I(list(smckillvec))
)

saveRDS(inputs, file = 'inputs.rds')

# Initialize list to store parasitaemia trajectories per child
parasitemia_storage <- vector("list", length = cohortparams$trial_ts + cohortparams$burnin)
N = cohortparams$N
trial_ts = cohortparams$trial_ts
burnin = cohortparams$burnin

susceptibles <- rep(TRUE, N)

for (t in 1:(trial_ts + burnin)){
  
  # Determine number of new infectious bites on time t
  n_infectious_bites <- rpois(1, lambda = N * p_bite_generic$prob_infectious_bite[t]) # or could do rbinom(1, N, p_bite)
  
  if (n_infectious_bites == 0) next
  
  # Sample random children to be bitten 
  bites <- sample(1:N, size = n_infectious_bites, prob = weights, replace = TRUE) #prob = weights,
  message('There are ', length(bites), " mosquito bites on time ", t, ' out of ', (trial_ts+burnin))
  
  bit_kids <- unique(bites)
  message('There are ', length(bit_kids), " unique kids bitten on time ", t, ' out of ', (trial_ts+burnin))
  
  # Update susceptibility vector so that recovered kids would be susceptible again
  if (exists("infection_records") && nrow(infection_records) > 0 & !allow_superinfections) {
    
    # Filter to only records with actual detections
    detected_records <- infection_records[!is.na(infection_records$detection_day), ]
    
    if (nrow(detected_records) > 0) {
      # Get the most recent detection for each person
      most_recent <- detected_records %>%
        group_by(rid) %>%
        slice_max(detection_day, n = 1) %>%  # Get most recent detection_day for each rid
        ungroup()
      
      # Find who should recover today (day t)
      recovery_today <- (most_recent$recovery_day + burnin) == t # add burnin because the recovery day calculation is the external time but for this siulation we want ot know if on time t they could be infected
      # if the recovery day equals today then these are the ones that are newly susceptible 
      if (any(recovery_today)) {
        recovering_kids <- most_recent$rid[recovery_today]
        susceptibles[recovering_kids] <- TRUE
      }
      
    }
  }
  
  # Avoid superinfection -- if a kid is bitten and they aren't susceptible they are not added to the list of bit_kids
  bit_kids_t <- c()
  j <- 1
  for (kid in bit_kids){
    if (susceptibles[kid]){ 
      bit_kids_t[j] <- kid # Add this kid to the filtered list if the kid is susceptible
      
      j <- j + 1 # move to the next position in the filtered list 
    }
  }
  message(str_glue('{length(bit_kids_t)} / {length(bit_kids)} bit kids were susceptible'))
  bit_kids <- bit_kids_t
  
  if(length(bit_kids) > 0) {
    # Estimate antibody level for each bitten child 
    if(t < burnin){
      # No interventions during burnin
      PEV_vec <- rep(0, length(bit_kids))
      SMC_vec <- rep(0, length(bit_kids))
      # these both need a placeholder value (won't matter if PEV = 0 or SMC = 0)
      t_since_vax_vec <- rep(0, length(bit_kids))
      t_toboost1_vec <- rep(500, length(bit_kids))
      t_toboost2_vec <- rep(1000, length(bit_kids))
      
      SMC_kill_vec <- rep(list(rep(0, floor(trial_ts/2 + burnin))), length(bit_kids))
      SMC_timev <- rep(list(0:(length(SMC_kill_vec[[1]])-1)), length = length(bit_kids))
    } else {# Find time in days between vaccination and infection
      
      # Vectorized intervention calculations
      kid_metadata <- metadata_child[metadata_child$rid %in% bit_kids, ]
      
      # Create named vectors
      pev_lookup <- setNames(kid_metadata$PEV, kid_metadata$rid)
      smc_lookup <- setNames(kid_metadata$SMC, kid_metadata$rid)
      vax_day_lookup <- setNames(kid_metadata$vaccination_day, kid_metadata$rid)
      
      # Map to bit_kids order
      PEV_vec <- pev_lookup[as.character(bit_kids)]
      SMC_vec <- smc_lookup[as.character(bit_kids)]
      t_since_vax_vec <- (t - burnin) - vax_day_lookup[as.character(bit_kids)] # here, a - value of vaxdaylookup indicates vaccination before the follow-up  
      
      # SMC
      # Make named vector of SMC and map to bit_kids order 
      smc_lookup <- setNames(kid_metadata$SMC, kid_metadata$rid)
      SMC_vec <- smc_lookup[as.character(bit_kids)]
      
      # Get the kill rate vectors for SMC 
      # Find the row index for the target rid
      smc_kill_vec_lookup <- setNames(kid_metadata$smckillvec, kid_metadata$rid)
      SMC_kill_vec <- smc_kill_vec_lookup[as.character(bit_kids)]
      # subset the kill rate vector to be from this external time 
      # divide max time by two because the vector is every two days 
      SMC_kill_vec <- lapply(SMC_kill_vec, function(smcvec){
        subset <- smcvec[floor((t + burnin) / 2) :length(smcvec)]
        return(subset)
      })
      
      SMC_timev <- rep(list(0:(length(SMC_kill_vec[[1]])-1)), length = length(bit_kids))#rep(list(0:(trial_ts/2 + burnin)), length = length(bit_kids))
      
      # SMC_kill_vec = unlist(inputs$smckillvec)
      # SMC_kill_vec = rep(list(SMC_kill_vec[floor((t + burnin) / 2) :length(SMC_kill_vec)]), bit_kids)
      # SMC_timev <- rep(list(seq(0, (length(SMC_kill_vec)-1))), length = bit_kids)

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
      infection_start_day = rep(t, length(bit_kids)),
      # season_start_day = season_start_day,
      # season_length = season_length, 
      # smc_interval = smc_interval,
      rid = bit_kids,
      SMC_time = I(SMC_timev),
      SMC_kill_vec = I(SMC_kill_vec),
      t_toboost1 = 365,
      t_toboost2 = 730
    )
    # Find time length until end of follow-up to run infection sim  
    # tt_until_end_cohort <- seq(0, ((trial_ts + burnin) - t)/divide, by = tstep)#
    # message('time until end of cohort, run infection sim for: ', tt_until_end_cohort, ' timesteps')
 
    outputs <- pmap(params_df, 
                    function(PEV_on, 
                             SMC_on, 
                             t_inf, 
                             infection_start_day,
                             # season_start_day, 
                             # season_length, 
                             # smc_interval, 
                             SMC_time,
                             SMC_kill_vec,
                             rid,
                             t_toboost1,
                             t_toboost2) {
                      result <- run_process_model(
        PEV_on = PEV_on,
        SMC_on = SMC_on,
        t_inf = t_inf,
        VB = VB,
        tt = tt,
        # max_SMC_kill_rate = max_SMC_kill_rate,
        # SMC_decay = SMC_decay,
        # smc_timing = smc_timing,
        infection_start_day = infection_start_day, # external time that infection begins 
        SMC_time = unlist(SMC_time),
        SMC_kill_vec = unlist(SMC_kill_vec),
        tboost1 = t_toboost1, 
        tboost2 = t_toboost2
        # season_start_day = season_start_day, # day of season start relative to Jan 1 of external year 
        # season_length = season_length, # ~4 months (120 days)
        # smc_interval= smc_interval # how often are SMC rounds (days)
      )
                      
      result$trajectory$rid <- rid
      
      return(result)
    })
    
    # Vectorized data frame creation for infection records
    new_records <- data.frame(
      time_ext = t - burnin,                                                      # external cohort time (t is the cohort time, then we wnat to scale to be + if after burnin)
      rid = bit_kids,                                                        # number of bit kids
      infectious_bite_day = (t - burnin) - t_liverstage,                          # bitten on day t, then assuming the liver stage takes t_liverstage days
      BSinfection_day = (t - burnin),                                             # after liver stage, the BS begins #+ t_liverstage
      threshold_day = sapply(outputs, function(x) x$threshold_day),               # days since BS starts that threshold is reached
      detection_day = (t - burnin) + sapply(outputs, function(x) x$threshold_day),# day in cohort simulation that threshold is reached, threshold is the day since BS infection that reaches threshold #+ t_liverstage 
      t_toreach_threshold = sapply(outputs, function(x) x$threshold_day) + t_liverstage,   # time to reach threshold value / detection since the bite  
      vaccination_day = if(t < burnin) rep(vax_day, length(bit_kids)) else kid_metadata$vaccination_day,   # day of vaccination relative to the start of follow-up (day 0 external time)
      prob_bite = p_bite_generic$prob_infectious_bite[t],
      recovery_day = ((t - burnin) + sapply(outputs, function(x) x$threshold_day)) + 12 - t_liverstage # day that the child would be 'recovered' if we assume that a child is treated and has a period of prophylaxis for 12 days after detection day
    ) 
    
    infection_records <- rbind(infection_records, new_records)
    
    # Only update the susceptibility vector if we are not allowing superinfections 
    if(!allow_superinfections){
      # mark as non-susceptible only those who have a valid detection_day
      detected_kids <- bit_kids[!is.na(new_records$detection_day)]
      # print(detected_kids)
      if (length(detected_kids) > 0) {
        susceptibles[detected_kids] <- FALSE
      } 
    }
    
    # Vectorized parasitemia storage creation
    parasitemia_data <- map2(outputs, bit_kids, function(output, kid) {
      
      output$trajectory %>%
        mutate(
          day1_BSinfection = t - burnin,#+ t_liverstage 
          detection_day = t + (output$threshold_day) - burnin,# do not need to multiply threshold day by 2 here since already done in run_process_model #+ t_liverstage 
          time_ext = time*2 + (t - burnin),# + infection_start_day/2,# external time should be dependent on when infection was, relative to external time(inf_start_day aka t); time (model time) has already been multiplied by 2 and related to infection time  #time*2 + (t - 1) - burnin,#+ t_liverstage
          intervention = metadata_child[metadata_child$rid == kid, ]$intervention,
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
  group_by(rid, day1_BSinfection) %>%
  mutate(child_dayinf = paste0(rid, ", day ", day1_BSinfection),
  ) %>%
  mutate(det = ifelse(!is.na(detection_day), 1, 0)) #%>% # this is the day of the cohort simulation that the infection is found --

saveRDS(parasitemia_df, file = 'outputs/parasitemia.rds')


# Add in relevant metadata to infection records df
child_counts <- metadata_child %>%
  distinct(rid, intervention) %>%
  count(intervention, name = "children_in_group")

infection_records2 <- infection_records %>% 
  filter(detection_day > 0 | is.na(detection_day)) %>% # this is the day of the cohort simulation that the infection is found --
  left_join(metadata_child %>% select(-c(vaccination_day)), by = c('rid')) %>%
  left_join(child_counts, by = 'intervention') %>%
  mutate(detectable = ifelse(is.na(threshold_day), 0, 1)) %>%
  group_by(intervention) %>%
  arrange(detection_day) %>%
  mutate(cumul_inf = cumsum(detectable))

saveRDS(infection_records2, file = 'outputs/infection_records.rds')


# Make plots 

# Parasitemia over time by child in sample of 100 random children per group
parasitemia_random <- ggplot(parasitemia_df[parasitemia_df$rid %in% (parasitemia_df %>% distinct(rid, intervention) %>% group_by(intervention) %>% slice_sample(n=50) %>% pull(rid)), ]) +
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
  geom_line(aes(x = time_ext/365.25, y = log10(numkillSMC), color = as.factor(detection_day))) +
  geom_line(aes(x = time_ext/365.25, y = prob_smckill, color = as.factor(detection_day))) +
  # geom_line(aes(x = time_ext/365.25, y = log10(parasites)/2), color = 'black') +
  geom_line(aes(x = time_ext/365.25, y = growth/2), color = 'pink', alpha = 0.5) +
  # geom_line(data = as.data.frame(rainfall), aes(x = t/365.25, y = profile*100)) +
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
