# RTS,S over time 

# parameters from  White 2015 Table S3

# Fitted parameter values: from White et al (2015) Lancet ID
ab_mu <- 621  # median of the geometric means of observed antibody titres
ab_mu_boost <- 277 # median of the boosted antibody titre
t_boost <- 548 #time of booster dose in days after third vaccine dose
ab_sigma <- 0.35 # observational variance of antibody titre (log-normal)
ab_sigma_boost <- 0.35
d1_mu <- 45 # half-life of the short-lived component of antibody response 
d1_sigma <- 16 # standard deviation of HL of short-lived response
d2_mu <- 591 # half-life of long-lived component of antibody response
d2_sigma <- 245 # standard deviation of HL of long-lived response
rho_mu <- 2.378 # proportion of short-lived component following primary schedule (logit-normal distribution)
rho_sigma <- 1.008 
rho_mu_boost <- 1.034 # proportion of short-lived component following booster (logit-normal distribution)
rho_sigma_boost <- 1.027

# FUnctions from Bob Verity (through https://github.com/mrc-ide/rtss_vacc_antibody_model/tree/main)
# draw from logitnormal distribution. 
# Input arguments include the mean and standard deviation of the RAW normal variate that will be converted to the logit scale. 
# Therefore the mean and standard deviation of the final random variable will not equal these values.
rlogitnorm <- function(n, mean_raw, sd_raw) {
  x <- exp(rnorm(n, mean_raw, sd_raw))
  x/(1+x)
}

# draw from lognormal distribution. 
# Input arguments include the desired mean and standard deviation of the final lognormal random variable; 
# the mean and standard deviation of the raw normal variate are calculated from these values.
rlnorm2 <- function(n, mean, sd) {
  meanlog <- log(mean^2/sqrt(mean^2+sd^2)) # mean of the lognormal distribution
  sdlog <- sqrt(log(1+(sd/mean)^2)) # standard deviation of the lognormal distribution
  rlnorm(n,meanlog,sdlog) # draw from the lognormal distribution
}

#' Antibody model originally from Hogan's implementation  https://github.com/mrc-ide/rtss_vacc_antibody_model/tree/main
#' modified to add second booster (assuming same parameters as first)
# antibody_titre <- function(t,
#                            ds_mu = 45,# half-life of the short-lived component of antibody response
#                            ds_sigma = 16,# standard deviation of HL of short-lived response
#                            dl_mu = 591,# half-life of long-lived component of antibody response
#                            dl_sigma = 245,# standard deviation of HL of long-lived response
#                            rho_mu = c(2.37832,1.034,1.034), # proportion of short-lived component following primary and booster schedule (logit-normal distribution)
#                            rho_sigma = c(1.00813,1.027,1.027), # sd of proportion of short-lived component following primary and booster schedule (logit-normal distribution)
#                            ab_mu = c(621, 277, 277), # median of geometric means of observed ab titres for primary series and booster ()
#                            ab_sigma = c(0.35, 0.35, 0.35), # observational variance of ab titre (log-normal)
#                            tboost1 = 364,
#                            tboost2 = 729
#                            ){
# 
#   # draw parameters from distributions
#   rho_prim <- rlogitnorm(1, mean_raw=rho_mu[1], sd_raw=rho_sigma[1])
#   rho_boost1 <- rlogitnorm(1, mean_raw=rho_mu[2], sd_raw=rho_sigma[2])
#   rho_boost2 <- rlogitnorm(1, mean_raw=rho_mu[3], sd_raw=rho_sigma[3])
#   rho <- c(rho_prim, rho_boost1, rho_boost2)
#   ds <- rlnorm2(1, mean=ds_mu, sd=ds_sigma) # half-life of short-lived component of antibody response
#   dl <- rlnorm2(1, mean=dl_mu, sd=dl_sigma) # half-life of long-lived component of antibody response
#   ab0 <- exp(rnorm(1, log(ab_mu[1])-ab_sigma[1]^2/2, sd=ab_sigma[1]))
#   ab0_boost1 <- exp(rnorm(1, log(ab_mu[2])-ab_sigma[2]^2/2, sd=ab_sigma[2]))
#   ab0_boost2 <- exp(rnorm(1, log(ab_mu[3])-ab_sigma[3]^2/2, sd=ab_sigma[3]))
#   csp_peak <- c(ab0, ab0_boost1, ab0_boost2)
# 
#   rs <- log(2) / ds
#   rl <- log(2) / dl
# 
#   t_toboost <- c(0, tboost1, tboost2)
# 
#   csp <- rep(NA, length(t))
#   for(i in seq_along(t_toboost)){
#     t_cur <- t - t_toboost[i]
#     index <- t >= t_toboost[i]
# 
#     csp[index] = csp_peak[i] * (rho[i] * exp(-rs * t_cur[index]) + (1 - rho[i]) * exp(-rl * t_cur[index]))
#   }
#   return(csp)
# }


#' Antibody model 
#' 
#' Following White et al (2015). Defaults are from White Table 3 and Table 1
#' 
#' @param t timesteps
#' @param csp_peak median GMTs of primary series and boosters, respectively 
#' @param rho proportion of short lived component following primary series and booster, respectively 
#' @param ds half life of short lived component in days
#' @param dl half life of long lived component in days
#' @param tboost timestep of booster dose
# antibody_titre <- function(t, 
#                            csp_peak = c(621, 277), 
#                            rho = c(0.88, 0.7), 
#                            ds = 45,   
#                            dl = 591,   
#                            tboost = 548){ 
#   rs <- log(2) / ds
#   rl <- log(2) / dl
#   
#   t_toboost <- c(0, tboost)
#   
#   csp <- rep(NA, length(t))
#   for(i in seq_along(t_toboost)){
#     t_cur <- t - t_toboost[i]
#     index <- t >= t_toboost[i]
#     
#     csp[index] = csp_peak[i] * (rho[i] * exp(-rs * t_cur[index]) + (1 - rho[i]) * exp(-rl * t_cur[index]))
#   }
#   return(csp)
# }  

# csp <- antibody_titre(t =c(1:(365*3)))
# csp <- antibody_titre_det(t =c(1:(365*3)))
# plot(csp)
#' Vaccine efficacy model https://github.com/mrc-ide/r21_vacc_antibody_model/blob/main/efficacy_model/R/vaccine_efficacy_model.R 
#' 
#' Following White et al (2015) Defaults are from White Table 3 (5-17 months)
#'
#' @param cspt Antibody titre
#' @param v_max Maximum efficacy against infection
#' @param alpha Shape parameter of dose–response curve
#' @param beta Scale parameter of dose–response curve (EU/mL)
#' these are all for 5-17 m category
#'
#' @return Vector of vaccine efficacy.
vaccine_eff <- function(csp, 
                        vmax = 0.93, 
                        alpha = 0.74,  
                        beta = 99.2){ 
  
  v_t = vmax * (1 - (1 / (1 + (csp / beta) ^ alpha)))
  return(v_t)
  
}

#' Vaccine efficacy per sporozoite following White et al.  2013
p_spz_surv <- function(ab, beta_ab = 6.62,
                       alpha_ab = 1.32){
  
  1 / (1 + (ab / beta_ab)^alpha_ab) 
  
}

# eff <- vaccine_eff(csp = csp)

# ggplot() +
#   geom_line(aes(x = 0:(365*3), y = eff))
# ggplot() +
#   geom_line(aes(x = 0:(365*3), y = csp)) + 
#   scale_y_log10()

# function modified from Nora Schmit https://github.com/mrc-ide/r21_vacc_antibody_model/blob/main/ab_model/R/antibody_model.R to include uncertainty
antibody_titre <- function(t, 
                           phase, 
                           peak1, peak2, peak3, 
                           duration1, duration2, 
                           rho1, rho2, rho3, 
                           t_boost1 = 364, t_boost2 = 729){
  # rho3 is same as rho2 unless we get more information that it is different 
  # duration1 and duration2 are mean, sd for short- and long-lived components 
  ds_draw = rlnorm2(1, mean=duration1[1], sd=duration1[2])
  dl_draw = rlnorm2(1, mean=duration2[1], sd=duration2[2])
  
  r1 <- log(2) / ds_draw
  r2 <- log(2) / dl_draw
  
  t[phase == 2] <-  t[phase == 2] - t_boost1
  t[phase == 3] <-  t[phase == 3] - t_boost2
  
  # rho1 and rho2 are mean, sd
  rho1_draw <- rlogitnorm(1, mean_raw=rho1[1], sd_raw=rho1[2])
  rho2_draw <- rlogitnorm(1, mean_raw=rho2[1], sd_raw=rho2[2])
  rho3_draw <- rlogitnorm(1, mean_raw=rho3[1], sd_raw=rho3[2])
  rho <- rep(rho1_draw, length(t))
  rho[phase == 2] <- rho2_draw
  rho[phase == 3] <- rho3_draw
  
  # draw peak csp
  # peak1, peak2, peak3 are edian of geometric means of observed ab titres and variance for primary, booster 1, booster 3
  peak1_draw <- exp(rnorm(1, log(peak1[1])-peak1[2]^2/2, sd=peak1[2]))
  peak2_draw <- exp(rnorm(1, log(peak2[1])-peak2[2]^2/2, sd=peak2[2]))
  peak3_draw <- exp(rnorm(1, log(peak3[1])-peak3[2]^2/2, sd=peak3[2]))
  peak <- rep(peak1_draw, length(t))
  peak[phase == 2] <- peak2_draw
  peak[phase == 3] <- peak3_draw
  
  ab <- peak * ((rho * exp(-r1 * t)) + ((1 - rho) *  exp(-r2 * t)))
  
  # Check for down-boosts
  last_phase_1 <- tail(ab[phase == 1], 1)
  if(last_phase_1 > peak2[1]){
    peak[phase == 2] <- last_phase_1
  }
  last_phase_2 <- tail(ab[phase == 2], 1)
  if(last_phase_2 > peak3[1]){
    peak[phase == 3] <- last_phase_2
  }
  ab <- peak * ((rho * exp(-r1 * t)) + ((1 - rho) *  exp(-r2 * t)))
  
  return(ab)
}

