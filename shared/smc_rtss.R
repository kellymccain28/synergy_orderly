# odin2 model for infection with SMC and RTSS 
# This model can be run either deterministically or stochastically 

## Core equations for transitions between compartments:
update(PB) <- if(PB_next < 1e-5) 0 else (PB_next)# - n_killSMC) * growth # first see how many die then multiply by growth
update(sc) <- 1 / (1 + (PB_next / pc)^kappac) # gets smaller as parasite density increases (could maybe set a minimum value?)
update(sm) <- if(time > 3) beta + (1 - beta) / (1 + (PBsum / pm)^kappam) + beta else 1
update(growth) <- m * sc * sm ## growth rate is modified by natural and general-adaptive immunity

## Initial states:
initial(PB) <- (mero_init / VB)
initial(sc) <- 1 # this means that multiplying by the growth rate keeps the growth rate as it would be without immunity
initial(growth) <- m_init
initial(sm) <- 1

# make variable to use what the value of PB will be at current timestep
PB_next <- (PB - n_killSMC) * growth #PB + n_grow - n_killSMC ##

# make variable to use what the value of PB will be at current timestep
# growthrate <- 1 - exp(-growth) * dt
# n_grow <- Poisson(((PB * VB) - n_killSMC) * growth) / VB

# Draw growth rate for each timestep (2 days) -- parameters from Challenger et al.
corr_m <- 0.63
mean_m <- 16
sd_m <- 8.7
min_m <- 1
max_m <- 35
mean_new <- mean_m + corr_m * (m - mean_m)
sd_new <- sd_m * sqrt(1 - corr_m^2)
m_init <- TruncatedNormal(mean = mean_m, sd = sd_m, min = min_m, max = max_m)
initial(m) <- m_init
update(m) <- TruncatedNormal(mean = mean_new, sd = sd_new, min = min_m, max = max_m)

# Make cumulative general adaptive immunity PC
update(PC) <- PC + if(PB_next > C) C else PB_next 
initial(PC) <- 0

# Make variable to keep track of  sum of parasitemia over timesteps t = 1 to t = t-4+1
n_dim <- 3
dim(buffer) <- n_dim
initial(buffer[]) <- 0

## Most recent data first
update(buffer[1]) <- if(PB_next > C) C else PB_next #if(PB_next > C) C else PB_next#
update(buffer[2:n_dim]) <- buffer[i - 1]

# Make total variable
initial(PBsum) <- 0
update(PBsum) <- if(time > 3) PC - sum(buffer) else 0


# SMC 
# Calculate day of year in cohort time to determine SMC curve
external_day <- if (SMC_on == 1) infection_start_day + time * 2 else 0 #infection_start_day/2 + time # this is the day of the infection in external time (not infection time)
day_of_year <- if (SMC_on == 1) ((external_day - 1) %% 365) + 1 else 0
season_end_day <- if (SMC_on == 1) season_start_day + season_length - 1 else 0 # season length is in days, but the timesteps are 2 days long
# Check if in season
in_season <- if (season_end_day <= 365) (
  if (day_of_year >= season_start_day && day_of_year <= season_end_day) 1 else 0
) else (
  if (day_of_year >= season_start_day || day_of_year <= (season_end_day - 365)) 1 else 0
)
initial(in_season_out) <- 0
update(in_season_out) <- in_season
# Days since season started (accounting for year wrap)
days_in_season_raw <- if (season_end_day <= 365) (
  day_of_year - season_start_day + 1
) else (
  if (day_of_year >= season_start_day) (
    day_of_year - season_start_day + 1
  ) else (
    day_of_year + (365 - season_start_day + 1)
  )
)
days_in_season <- max(1, days_in_season_raw)
# Time since last SMC dose
smc_dose_day <- floor((days_in_season - 1) / (smc_interval))
time_since_smc <- (days_in_season - 1) - smc_dose_day * (smc_interval)

# SMC kill rate
SMC_kill_rate <- if (SMC_on == 1) in_season * max_SMC_kill_rate * exp(-(time_since_smc / lambda)^kappa) else 0

# probability of parasite being killed by SMC
lambda <- parameter(17) #(scale)
kappa <- parameter(0.28) #(shape)
# SMC_kill_rate <- if (time >= smc_timing) max_SMC_kill_rate * exp(-((time - smc_timing) / lambda)^kappa) else 0 # hill
# SMC_kill_rate <- if (time >= smc_timing) max_SMC_kill_rate * exp(-SMC_decay * (time - smc_timing)) else 0#max_SMC_kill_rate * exp(-SMC_decay * dt) # ## time-varying SMC kill rate (exponential decay) # exp
p_killSMC <-  if (SMC_on == 1) 1 - exp(-SMC_kill_rate * dt) else 0
n_killSMC <- if (SMC_on == 1) Binomial(PB * VB, p_killSMC) / VB else 0

# print("external_day: {external_day}, day_of_year: {day_of_year}, in_season: {in_season}, days_in_season: {days_in_season}, smc_dose_day: {smc_dose_day}, time_since_smc: {time_since_smc}, SMC_kill_rate: {SMC_kill_rate}")

initial(SMC_kill_rateout) <- 0
update(SMC_kill_rateout) <- SMC_kill_rate
initial(prob_smckill) <- 0
update(prob_smckill) <- p_killSMC
initial(numkillSMC) <- 0
update(numkillSMC) <- n_killSMC


######## Immunity 10.1038/s41467-017-01352-3

# Draw from log-normal distribution for innate response 
kappac <- 3
kc <- 0.164
mu_pc <- 4.79 * log(10) # mean of LN distribution
sig_pc <- 1.2 # sd of log normal distribution
rand_pc <- exp(TruncatedNormal(mean = mu_pc,
                               sd = sig_pc,
                               min = 0,
                               max = log(10^5.5)))
pc <- kc * rand_pc

# Draw from Gompertz distribution for general adaptive response 
kappam <- 1                 #
C <- 1                      # regulates growth rate of general adaptive immunity
beta <- 0.01                # controls max efficacy of the general adaptive immunity
km <- 0.021                 #
alpha_pm <- 0.0311 #1       # first parm for Gompertz distribution to control gen adaptive response
theta_pm <- 0.0004 #2       # second " "

unif <- Uniform(0,1)
pm_rand <- log(1 - (1 / theta_pm) * log(1 - unif)) * (1 / alpha_pm)
pm <- km * pm_rand

########

# Draw number of surviving sporozoites and subsequent merozoites initiating blood stage infection
#Dose response
ab <- if(PEV_on == 1) ab_user else 0
DR <- 1 / (1 + (ab / beta_ab)^alpha_ab) # prob of survival of single spz

# Parameters for spz model initial merozoites 
# estimated and fixed parameters  
n <- 150/5 #n, mean number of successful spz per challenge / 5 bites; neg bin
sigma_n <- 194/5 #, sigman sd of number of successful spz per challenge / 5 bites
mu <- 2136 #30000 #10.1371/journal.pcbi.1005255 as assumed by Hayley #2136 # mean number of merozoites released per sporozoite in Michael's model; gamma distributed
sigma_mu <- 4460 #71427 #4460 from Michael's model # sd of number of merozoites released per sporozoite; gamma distributed
beta_ab <- 6.62 #8639 # from hayley's #6.62 # anti-CSP titre for 50% reduction in spz survival prob microgram/mL
alpha_ab <- 1.32 #1.53 # from hayley's #1.32 # shape parameter for antibody dose-response

# Parameters for Negative Binomial distribution
r <- n^2 / (sigma_n^2 - n)
p <- r / (n*DR + r) 

# Draw number of successful sporozoites
kspz <- NegativeBinomial(r, p)

# Parameters for Gamma distribution
theta <- sigma_mu^2 / mu  # Scale parameter (theta)
gamma_shape <- kspz * mu / theta
# find total merozoites for from k sporozoites
mero_init <- Gamma(shape = gamma_shape, scale = theta) 

## User defined parameters - default in parentheses:
# m <- parameter(20) # average PMR over two days or 48 hours (initial value)
PEV_on <- parameter()
SMC_on <- parameter()
ab_user <- parameter()
VB <- parameter()
max_SMC_kill_rate <- parameter(1.5)
# SMC_decay <- parameter(0.05) # decay of SMC kill rate
# smc_timing <- parameter(0) # timing of receipt of last dose of SMC relative to start of simulation

# SMC parameters
infection_start_day <- parameter(0) # external time that the infection began -- default is 0 (if running outside of cohort )
season_start_day <- parameter(10)    # Day 180 = ~July 1st
season_length <- parameter(120)       # 4 months (120 days)
smc_interval <- parameter(30)         # SMC every 30 days

