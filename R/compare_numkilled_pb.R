# test that num killed in smc is lower than in timestep before

smc_dose_days <- 0
smc_dose_days <- c(seq(season_start_day, season_start_day + 120 - 1, 30),
                   seq(season_start_day + 366, season_start_day + 366 + 120 - 1, 30),
                   seq(season_start_day + 365*2, season_start_day  + 366*2 + 120 - 1, 30))
# will output kill vector in length ts
tssmc = 365*3/2
smc_killvec <- unlist(get_smc_vectors(smc_dose_days = smc_dose_days,
                                      ts = tssmc, 
                                      max_SMC_kill_rate = runpars$max_SMC_kill_rate,
                                      lambda = runpars$lambda, 
                                      kappa = runpars$kappa))
smckillvec_subset <- smc_killvec[floor((runpars$inf_start) / 2) :length(smc_killvec)]
smckillvec_subset <- list(c(smckillvec_subset, rep(0, tssmc-length(smckillvec_subset))))
smckilltime <- seq(0, length(smckillvec_subset[[1]])-1,1)
# Without vaccination but with SMC ----
smc <- run_model(n_particles = runpars$n_particles,
                 n_threads = 4L,
                 PEV_on = 0,
                 SMC_on = 1,
                 tt= tt,
                 det_mode = FALSE,
                 t_inf_vax = runpars$t_inf_vax,
                 VB = VB,
                 SMC_time = smckilltime,
                 SMC_kill_vec = smckillvec_subset,
                 infection_start_day = runpars$inf_start
) 

output <- data.frame(timestep = 1:200, 
                     pb = smc$PB[1,], 
                     growth = smc$growth[1,], 
                     prob_smckill = smc$prob_smckill[1,], 
                     numkillSMC = smc$numkillSMC[1,], 
                     pb_for_bin = smc$PB_for_binomial[1,], 
                     bin_size = smc$binomial_size[1,])
output <- output %>% 
  mutate(killedlessthanpb = ifelse(numkillSMC <= pb_for_bin, 1, 0 ),
         laggedpb = lag(pb),
         killedlessthanpblag = ifelse(numkillSMC <= laggedpb, 1, 0 ))
write.csv(output, 'output_testing.csv')


smc_dose_days <- 0
smc_dose_days <- c(seq(season_start_day, season_start_day + 120 - 1, 30),
                   seq(season_start_day + 366, season_start_day + 366 + 120 - 1, 30),
                   seq(season_start_day + 365*2, season_start_day  + 366*2 + 120 - 1, 30))
# will output kill vector in length ts
tssmc = 365*3/2
smc_killvec <- unlist(get_smc_vectors(smc_dose_days = smc_dose_days,
                                      ts = tssmc, 
                                      max_SMC_kill_rate = runpars$max_SMC_kill_rate,
                                      lambda = runpars$lambda, 
                                      kappa = runpars$kappa))
smckillvec_subset <- smc_killvec[floor((runpars$inf_start) / 2) :length(smc_killvec)]
smckillvec_subset <- list(c(smckillvec_subset, rep(0, tssmc-length(smckillvec_subset))))
smckilltime <- seq(0, length(smckillvec_subset[[1]])-1,1)
# Without vaccination but with SMC ----
smc <- run_model(n_particles = runpars$n_particles,
                 n_threads = 4L,
                 PEV_on = 0,
                 SMC_on = 1,
                 tt= tt,
                 det_mode = FALSE,
                 t_inf_vax = runpars$t_inf_vax,
                 VB = VB,
                 SMC_time = smckilltime,
                 SMC_kill_vec = smckillvec_subset,
                 infection_start_day = runpars$inf_start
) %>%
  format_data(tt= tt,
              infection_start_day = runpars$inf_start,
              n_particles = runpars$n_particles) %>%
  make_plots()
dfsmc <- smc[[6]]
dfsmc <- dfsmc %>%
  group_by(run) %>%
  arrange(time_withinhost) %>%
  mutate(#killedlessthanpb = ifelse(nkillsmc <= pb_for_bin, 1, 0 ),
         laggedpb = lag(parasites),
         killedlessthanpblag = ifelse(nkillsmc <= laggedpb, 1, 0 ))

table(dfsmc$killedlessthanpblag)

ggplot(dfsmc) + 
  geom_line(aes(x  = time, y = smc_prob)) +
  # geom_line(aes(x = time, y = nkillsmc/10000, color = run)) +
  # geom_line(aes(x = time+2, y = parasites/10000, color = run), linetype = 2) + 
  geom_line(aes(x = time, y = nkillsmc/laggedpb, color = run)) +
  scale_y_log10()
