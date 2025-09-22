# mini version of SMC delivery with vector of smc timings 

mod1 <- odin2::odin({
  
  SMC_time = parameter(constant= TRUE)
  SMC_kill_vec = parameter(constant = TRUE)
  dim(SMC_time, SMC_kill_vec) <- parameter(rank = 1)
  
  SMC_kill_rate = interpolate(SMC_time, SMC_kill_vec, 'constant')
  
  p_killSMC <-  1 - exp(-SMC_kill_rate * dt)
  
  initial(SMC_kill_rateout) <- 0
  update(SMC_kill_rateout) <- SMC_kill_rate
  initial(prob_smckill) <- 0
  update(prob_smckill) <- p_killSMC
  
})

max_SMC_kill_rate = 18
lambda =17 
kappa = 0.28
smc_timing <- c(5,15,30,55)
tt = seq(0, 100,1)
time_since_smc <- sapply(tt, function(d) {
  last_dose <- max(smc_timing[smc_timing <= d])
  ifelse(is.finite(last_dose), d - last_dose, NA)
})
# time_since_smc <- lapply(smc_timing, calc_time_since_dose, days = tt)

smckillvec = max_SMC_kill_rate * exp(-(time_since_smc / lambda)^kappa) 
pars <- list(
  SMC_kill_vec = smckillvec,
  SMC_time = tt
)
sys <- dust_system_create(mod1, 
                          pars, 
                          n_particles = 1,
                          n_threads = n_threads,
                          # seed = 12345L,
                          deterministic = FALSE)
# Set initial values using initial() equations in model
dust_system_set_state_initial(sys)

# Run model 
# tt <- seq(1,100,1)
bs_model <- dust_system_simulate(sys, tt)

out <- dust_unpack_state(sys, bs_model)
out

plot(out$SMC_kill_rateout)
plot(out$prob_smckill)
