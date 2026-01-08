# odin2 model for infection with SMC and RTSS 
# This model can be run either deterministically or stochastically 

## Core equations for transitions between compartments:
update(PB) <- if(PB_next < 1e-5) 0 else (PB_next)# - n_killSMC) * growth # first see how many die then multiply by growth
update(sc) <- 1 / (1 + (PB_next / pc)^kappac) # gets smaller as parasite density increases (could maybe set a minimum value?)
update(sm) <- if(time > 3) (1 - beta) / (1 + (PBsum / pm)^kappam) + beta else 1
update(sv) <- if(time > 3) 1 / (1 + (PBvsum / pv) ^ kappav) else 1
update(growth) <- m * sc * sm * sv## growth rate is modified by natural and general-adaptive immunity and var-specific immunity

## Initial states:
initial(PB) <- (mero_init / VB)
initial(sc) <- 1 # this means that multiplying by the growth rate keeps the growth rate as it would be without immunity
initial(growth) <- m_init
initial(sm) <- 1
initial(sv) <- 1

# make variable to use what the value of PB will be at current timestep
PB_next <- (PB - n_killSMC) * growth #PB + n_grow - n_killSMC ##

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

# Make total variable -- cumulative minus the parasitemia from the last few timesteps
initial(PBsum) <- 0
update(PBsum) <- if(time > 3) PC - sum(buffer) else 0
# update(PBsum) <- if(time > 3) sum(buffer) else 0



# Make total paraistemia variable for var-specific immunity Pv, bounds checking and sum calculation
# No sum before t=4, # Invalid range if ft>tminus4,  # Valid range - sum from f_t_timestep to t_minus_4_timestep
update(sv_timestep_range) <- t_minus_4 - f_t
initial(sv_timestep_range) <- 0
update(PBvsum) <- if (time <= 4) 0 else if (sv_timestep_range <= 0) 0 else sum(p_history[t_minus_4_ts:f_t_ts])
initial(PBvsum) <- 0

# Make array variable to keep track of paraistemia over gradually expanding time window 
# lower timestep is 4 timesteps back in the arrray that is shifting over time 
update(t_minus_4) <- if (time <= 4) 1 else time - 4
initial(t_minus_4) <- 1
update(t_minus_4_ts) <- time - t_minus_4 + 1 # get the index of the array p_history for lowe rrange
initial(t_minus_4_ts) <- 1
# upper timestep 
update(f_t) <- if (time <= 4) 1 else max(1, ceiling(lambdav^(time - 4)* (time - 4)))
initial(f_t) <- 1
update(f_t_ts) <- time - f_t + 1# get the index of the array p_history for upper range
initial(f_t_ts) <- 1

n_pv_dim <- tt # this is number of timesteps for the within-host model
dim(p_history) <- n_pv_dim
initial(p_history[]) <- 0
# more recent data are first
update(p_history[1]) <- PB_next
update(p_history[2:n_pv_dim]) <- p_history[i - 1]# for current timestep, add in the current parasitemia


# SMC 
# # If we want to run the model in generic SMC mode (not specific to individuals): 
# # Calculate day of year in cohort time to determine SMC curve
# external_day <- if (SMC_on == 1) infection_start_day + time * 2 else 0 #infection_start_day/2 + time # this is the day of the infection in external time (not infection time)
# day_of_year <- if (SMC_on == 1) ((external_day - 1) %% 365) + 1 else 0
# season_end_day <- if (SMC_on == 1) season_start_day + season_length - 1 else 0 # season length is in days, but the timesteps are 2 days long
# # Check if in season
# in_season <- if (season_end_day <= 365) (
#   if (day_of_year >= season_start_day && day_of_year <= season_end_day) 1 else 0
# ) else (
#   if (day_of_year >= season_start_day || day_of_year <= (season_end_day - 365)) 1 else 0
# )
# initial(in_season_out) <- 0
# update(in_season_out) <- in_season
# # Days since season started (accounting for year wrap)
# days_in_season_raw <- if (season_end_day <= 365) (
#   day_of_year - season_start_day + 1
# ) else (
#   if (day_of_year >= season_start_day) (
#     day_of_year - season_start_day + 1
#   ) else (
#     day_of_year + (365 - season_start_day + 1)
#   )
# )
# days_in_season <- max(1, days_in_season_raw)
# # Time since last SMC dose
# smc_dose_day <- floor((days_in_season - 1) / (smc_interval))
# time_since_smc <- (days_in_season - 1) - smc_dose_day * (smc_interval)
# # SMC kill rate
# SMC_kill_rate <- if (SMC_on == 1) max_SMC_kill_rate * exp(-(time_since_smc / lambda)^kappa) else 0

# probability of parasite being killed by SMC
# lambda <- parameter(17) #(scale)
# kappa <- parameter(0.28) #(shape)
# SMC_kill_rate <- if (time >= smc_timing) max_SMC_kill_rate * exp(-((time - smc_timing) / lambda)^kappa) else 0 # hill
# SMC_kill_rate <- if (time >= smc_timing) max_SMC_kill_rate * exp(-SMC_decay * (time - smc_timing)) else 0#max_SMC_kill_rate * exp(-SMC_decay * dt) # ## time-varying SMC kill rate (exponential decay) # exp
p_killSMC <-  if(SMC_on == 1) 1 - exp(-SMC_kill_rate * dt) else 0
n_killSMC <- Binomial(PB * VB, p_killSMC) / VB 

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

# EVSR immune response 
kappav <- 1
pv <- 10100 # controls strength of EVSR
lambdav <- 0.9996 # part of f(t) that governs duration of EVSR memory
# ts = seq(1,100,1)
# plot(ceiling(lambdav^(ts-4)*(ts-4)))


########
# https://github.com/mrc-ide/rtss_vacc_antibody_model
# Draw number of surviving sporozoites and subsequent merozoites initiating blood stage infection
# Dose response
ab <- if(PEV_on == 1) ab_user else 0
DR <- (1 / (1 + (ab / beta_ab)^alpha_ab)) # prob of survival of single spz from page 4 of white 2013

# Parameters for spz model initial merozoites 
# estimated and fixed parameters  
n <- 150 #n, mean number of successful spz per challenge ; neg bin
sigma_n <- 194 #, sigman sd of number of successful spz per challenge 
mu <- 2136  # mean number of merozoites released per sporozoite in Michael's model; gamma distributed
sigma_mu <- 4460  #from Michael's model # sd of number of merozoites released per sporozoite; gamma distributed
beta_ab <- parameter(5.83)#6.62) # anti-CSP titre for 50% reduction in spz survival prob microgram/mL
alpha_ab <- parameter(1.38)#1.32)  # shape parameter for antibody dose-response
vmin <- parameter(0) # minimum survival probability  (addition to white model to reduce effectiveness of the )


# Parameters for Negative Binomial distribution 
# adapted from : https://github.com/ht1212/quality_quantity_modelling/blob/master/R3_Efficacy_Function_IR/3_VE_per_Sporozoite
r <- n^2 / (sigma_n^2 - n) 
p <- vmin + (1 - vmin) * n*DR / (n*DR + r) # p here is probability that a sporozoite dies (which is 1-prob of survival, so need to use 1-p in negbinom call) - this is from https://github.com/ht1212/quality_quantity_modelling/blob/master/R2_Model_Fitting/1_MCMC_Models 
# r / (n*DR + r) - this is 1- n*DR / (n*DR + r) https://stat.ethz.ch/R-manual/R-devel/library/stats/html/NegBinomial.html
# the negbin function returns the number of failures before r successes, so for it to output the number of successful spz we need to invert p (I think)

# Draw number of successful sporozoites
num_bites <- parameter(1) # default is a single bite 
kspz1 <- NegativeBinomial(r / 5, 1-p) #/ 5 bites
kspz2 <- if(num_bites > 1) NegativeBinomial(r / 5, 1-p) else 0 # if there is a second bite on the same day, use this 
kspz3 <- if(num_bites > 2) NegativeBinomial(r / 5, 1-p) else 0 # if 3 bites on teh same day 
kspz <- kspz1 + kspz2 + kspz3 # add up the number of spz for each bite

# Parameters for Gamma distribution
theta <- sigma_mu^2 / mu  # Scale parameter (theta)
gamma_shape <- kspz * mu / theta
# find total merozoites for from k sporozoites
mero_init <- Gamma(shape = gamma_shape, scale = theta) 
initial(mero_init_out) <- mero_init
update(mero_init_out) <- if(time == 0) mero_init else 0
## User defined parameters - default in parentheses:
# m <- parameter(20) # average PMR over two days or 48 hours (initial value)
PEV_on <- parameter()
SMC_on <- parameter()
ab_user <- parameter()
VB <- parameter()
tt <- parameter(100, constant = TRUE)

# SMC parameters
SMC_time <- parameter(constant = TRUE) # timesteps over which the SMC kill rate should be interpolated
SMC_kill_vec <- parameter(constant = TRUE) # kill rate vector for individual child 
dim(SMC_time, SMC_kill_vec) <- parameter(rank = 1) # this means that both of the above 2 parameters are vectors 

SMC_kill_rate <- interpolate(SMC_time, SMC_kill_vec, 'constant')
