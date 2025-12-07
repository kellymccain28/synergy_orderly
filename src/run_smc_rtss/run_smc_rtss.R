library(tidyverse)
library(ggplot2)
library(odin2)
library(dust2)
library(cowplot)
library(glue)
library(orderly)

# orderly_strict_mode()
runpars <- orderly_parameters(n_particles = NULL,
                   n_threads = NULL,
                   t_inf_vax = NULL,# if before infection, then + number; time of infectious bite relative to vaccination, influences AB titre and also influences the length of SMC kill vec -- t_inf/2:end
                   ts = NULL,
                   tstep = NULL,
                   max_SMC_kill_rate = 2.333333,
                   lambda = 16.66667, 
                   kappa = 0.2222222, 
                   # SMC_decay = NULL,
                   season_start_day = NULL, # influences when SMC is delivered relative to the start of the infection/sim 
                   # season_length = NULL,
                   # smc_interval = NULL,
                   inf_start = NULL)# This is used in the formatting to shift the timesteps since start of BS for individual runs, and to subset the smc vector 
t_inf_vax = runpars$t_inf_vax
ts = runpars$ts
tstep = runpars$tstep
n_particles = runpars$n_particles
n_threads = runpars$n_threads
max_SMC_kill_rate = runpars$max_SMC_kill_rate
lambda = runpars$lambda
kappa = runpars$kappa
season_start_day = runpars$season_start_day
inf_start = runpars$inf_start
VB = 1e6
tt <- seq(1, runpars$ts, by = runpars$tstep)#
threshold <- 5000

orderly_shared_resource("smc.R",
                        "rtss.R",
                        "smc_rtss.R",
                        "helper_functions.R")

gen_bs <- odin2::odin("smc_rtss.R")
source("rtss.R")
source("helper_functions.R")

set.seed(1234)

# Stochastic runs ----
# Without vaccination nor SMC ----
# n_particles = 500L
# inf_start = 0-- this is when the infection shouldl start relative to the start of follow-up time (more relevant for cohort), can change timing of SMC
nothing <- run_model(n_particles = n_particles,
                     n_threads = 4L,
                     PEV_on = 0,
                     SMC_on = 0,
                     tt= tt,
                     SMC_time = seq(0,100,1),
                     SMC_kill_vec = rep(0,101),
                     t_inf_vax = t_inf_vax,
                     infection_start_day = inf_start,
                     VB = VB,
                     det_mode = FALSE) %>%
  format_data(tt= tt,
              infection_start_day = inf_start,
              n_particles = n_particles) %>%
  make_plots() 
prbc <- nothing[[1]] + xlim(c(0,200))
innate <- nothing[[2]]
genad <- nothing[[3]]
varspec <- nothing[[4]]
growthr <- nothing[[5]]
# nothing[[6]]
mplot <- nothing[[7]]
nothing[[11]] 

plot_grid(prbc+labs(caption = NULL, x = 'Days'),
          growthr+labs(caption = NULL, x = 'Days'),
          innate+labs(caption = NULL, x = 'Days'),
          genad+labs(caption = NULL, x = 'Days'),
          varspec+labs(caption = NULL, x = 'Days'),
          mplot+labs(caption = NULL, x = 'Days')
) # in caption, days since start of blood-stage 
ggsave('R:/Kelly/synergy_orderly/figures/immunity_plot.pdf')
nothingplt <- plot_grid(nothing[[1]] +labs(caption = NULL, x = 'Days since start of blood stage'), 
                        nothing[[2]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                        nothing[[3]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                        nothing[[4]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                        nothing[[6]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                        nothing[[7]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                        nothing[[8]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                        nothing[[9]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                        nrow = 4) 
dfnothing <- nothing[[6]] %>%
  mutate(scen = 'none')
# % of runs without infection at different times
# dfnothing %>% filter(time == 1, parasites < 1e-5) %>% count() / n_particles
# dfnothing %>% filter(time == max(time), parasites < 10) %>% count() / n_particles # higher means more runs are not infected
# dfnothing %>% filter(time == 15, parasites < 10) %>% count() / n_particles

# With vaccination ----
vax <- run_model(n_particles = runpars$n_particles,
                 n_threads = 4L,
                 PEV_on = 1,
                 SMC_on = 0,
                 t_inf_vax = runpars$t_inf_vax,
                 infection_start_day = runpars$inf_start,
                 SMC_time = seq(0,100,1),
                 SMC_kill_vec = rep(0,101),
                 alpha_ab = 1.32, 
                 beta_ab = 6.62,
                 vmin = 0, # higher vmin means that the probability of a sporozoite surviving is higher, so should have more infections
                 tt= tt,
                 VB = VB,
                 det_mode = FALSE) %>%
  format_data(tt= tt,
              infection_start_day = runpars$inf_start,
              n_particles = runpars$n_particles) %>%
  make_plots()
vax[[1]]# + xlim(c(0,100))
vax[[2]]
vax[[3]]
vax[[4]]
vax[[7]]
vax[['meroinit']]
# vax[[9]]
vaxplt <- plot_grid(vax[[1]] +labs(caption = NULL, x = 'Days since start of blood stage'), 
                    vax[[2]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                    vax[[3]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                    vax[[4]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                    vax[[6]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                    vax[[7]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                    vax[[8]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                    vax[[9]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                    nrow = 4) 
dfvax <- vax[[6]] %>%
  mutate(scen = 'vax')
# % of runs with without at day 60 
# dfvax %>% filter(time == max(time), parasites < 10) %>% count() / n_particles
# dfvax %>% filter(time == 15, parasites < 10) %>% count() / n_particles

# Calculate SMC kill vector ----
# intervalsmc_days = 30 
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
smc[[1]] + xlim(c(0,365))# plotting time which is in 2 day timesteps
# smc[[2]] #+ xlim(c(0,100))
# smc[[3]]
# smc[[4]]
# smc[[6]]
# smc[[7]]# + labs(x = 'Years', caption = NULL) #+ xlim(c(0,1))
smc[[8]]+ xlim(c(0,400))
smc[[9]]#+ xlim(c(0,100)) + ylim(c(0, 1))
smc[[10]]#+xlim(c(0,100))
smcplt <- plot_grid(smc[[1]]+labs(caption = NULL, x = 'Days since start of blood stage'), 
                    smc[[2]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                    smc[[3]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                    smc[[4]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                    smc[[6]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                    smc[[7]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                    smc[[8]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                    smc[[9]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                    nrow = 4) 
dfsmc <- smc[[6]] %>%
  mutate(scen = 'smc')
smc[['meroinit']] + labs(caption = NULL)

colorsplot <- c("#4123E8","#A6BF19")
parasites <- ggplot(dfsmc) + 
  geom_line(aes(x = time_withinhost2, y = parasites, group = run, color = detectable), alpha = 0.6, linewidth = 0.6) + #, color = cleared
  geom_line(aes(x = time_withinhost2, y = median_parasites, color = detectable), linewidth = 0.6)+
  geom_hline(aes(yintercept = threshold), linetype = 2, color = 'darkred', linewidth = 1) + # this is the detection limit (followiung Challenger et al.)
  geom_hline(aes(yintercept = 1e-5), linetype = 2, color = 'darkgreen', linewidth = 1) + # this is the clearance threshold
  scale_y_log10(labels = scales::label_log(),
                guide = "axis_logticks") +
  scale_x_continuous(limits = c(0, max(dfsmc$time_withinhost2)))+
  scale_color_manual(values = colorsplot) +
  labs(x = 'Days since start of blood stage',
       y = 'Parasites/uL') + 
  theme_bw() + xlim(c(0,400))+
  theme(legend.position = 'none')
smckillplot <- ggplot(dfsmc) + 
  geom_line(aes(x = time_withinhost2, y = smcrate, group = run, color = detectable), alpha = 0.7, linewidth = 0.8) + #
  scale_x_continuous(limits = c(0, max(dfsmc$time_withinhost2)))+
  labs(x = 'Days since start of blood stage',
       y = 'SMC per-parasite clearance rate\nper 2-day timestep') + 
  scale_color_manual(values = c('grey','maroon') ) +
  theme_bw() + xlim(c(0,400)) +
  theme(legend.position = 'none') 
smcprobplot <- ggplot(dfsmc) + 
  geom_line(aes(x = time_withinhost2, y = smc_prob, group = run, color = detectable), alpha = 0.7, linewidth = 0.6) + #
  scale_x_continuous(breaks = c(0, 7, seq(14, max(dfsmc$time_withinhost2), 28)),
                     limits = c(0, max(dfsmc$time_withinhost2)))+
  labs(x = 'Days since start of blood stage',
       y = 'SMC per parasite/uL\nclearance probability') + 
  scale_color_manual(values = c('grey','darkorchid4') ) +
  theme_bw() + 
  theme(legend.position = 'none') + xlim(c(0,400)) 
numkillsmcplot <- ggplot(dfsmc ) + 
  geom_line(aes(x = time_withinhost2, y = nkillsmc, group = run, color = detectable), alpha = 0.7, linewidth = 0.6) + #
  scale_x_continuous(breaks = c(0, 7, seq(14, max(dfsmc$time_withinhost2), 14)),
                     limits = c(0, max(dfsmc$time_withinhost2)))+
  geom_hline(aes(yintercept = threshold), linetype = 2, color = 'darkred', linewidth = 1) + # this is the detection limit (followiung Challenger et al.)
  geom_hline(aes(yintercept = 1e-5), linetype = 2, color = 'darkgreen', linewidth = 1) + # this is the clearance threshold
  labs(x = 'Days since start of blood stage',
       y = 'Parasites/uL cleared by SMC') + 
  scale_color_manual(values = colorsplot) +
  theme_bw() + 
  scale_y_log10(labels = scales::label_log(),
                guide = "axis_logticks") +
  theme(legend.position = 'none') + xlim(c(0,400)) 
  
smcplot <- plot_grid(parasites, smckillplot, numkillsmcplot, smcprobplot,
          labels = 'AUTO')
saveRDS(smc[[6]], paste0("R:/Kelly/synergy_orderly/figures/data/smc_data", Sys.Date(), '.rds'))
ggsave('R:/Kelly/synergy_orderly/figures/smcdynamics_plot.pdf', height = 6, width = 8)
# % of runs with infection at day 60 
# dfsmc %>% filter(time == max(time), parasites < 10) %>% count() / n_particles
# dfsmc %>% filter(time == 15, parasites < 10) %>% count() / n_particles

# Both vaccination and SMC
vaxSMC <- run_model(n_particles = runpars$n_particles,
                    n_threads = 4L,
                    PEV_on = 1,
                    SMC_on = 1,
                    tt= tt,
                    det_mode = FALSE,
                    t_inf_vax = runpars$t_inf_vax,
                    VB = VB,
                    alpha_ab = 1.32,
                    beta_ab = 6.62,
                    # vmin = 0, # higher vmin means that the probability of a sporozoite surviving is higher, so should have more infections
                    SMC_time = smckilltime,
                    SMC_kill_vec = smckillvec_subset,
                    infection_start_day = runpars$inf_start
) %>%
  format_data(tt= tt,
              infection_start_day = runpars$inf_start,
              n_particles = runpars$n_particles) %>%
  make_plots() 
vaxSMC[[1]]# + xlim(c(0,100))
# vaxSMC[[2]]
# vaxSMC[[3]]
# vaxSMC[[4]]
# vaxSMC[[7]]
vaxSMC[[11]] + labs(caption = NULL)
vaxSMCplt <- plot_grid(vaxSMC[[1]] +labs(caption = NULL, x = 'Days since start of blood stage'), 
                       vaxSMC[[2]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                       vaxSMC[[3]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                       vaxSMC[[4]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                       vaxSMC[[6]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                       vaxSMC[[7]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                       vaxSMC[[8]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                       vaxSMC[[9]]+labs(caption = NULL, x = 'Days since start of blood stage'),
                       nrow = 4) 
dfvaxsmc <- vaxSMC[[6]] %>%
  mutate(scen = 'vaxsmc')
# % of runs with infection at day 60 
# dfvaxsmc %>% filter(time == max(time), parasites < 10) %>% count() / n_particles
# dfvaxsmc %>% filter(time == 15, parasites < 10) %>% count() / n_particles


ggsave(filename = "nointerventions_stoch.png", plot = nothingplt, width = 11, height = 16)
ggsave(filename = "vaccine_stoch.png", plot = vaxplt, width = 11, height = 16)
ggsave(filename = "smc_stoch.png", plot = smcplt, width = 11, height = 16)
ggsave(filename = "vaxsmc_stoch.png", plot = vaxSMCplt, width = 11, height = 16)

# saveRDS(dfnothing, file = str_glue("df_nothing_", runpars$n_particles,".rds"))
# saveRDS(dfvax, file = str_glue("df_vax_", runpars$n_particles,".rds"))
# saveRDS(dfsmc, file = str_glue("df_smc_", runpars$n_particles,".rds"))
# saveRDS(dfvaxsmc, file = str_glue("df_vaxsmc_", runpars$n_particles,".rds"))

# poster
n <- nothing[[1]] + xlim(c(0,400)) + scale_color_manual(values = c("#F8766D")) +labs(caption = NULL) + theme_bw(base_size = 16) + theme(legend.position = 'none')
v <- vax[[1]] + xlim(c(0,400)) + scale_color_manual(values = c('#7CAE00'))+labs(caption = NULL)+ theme_bw(base_size = 16) + theme(legend.position = 'none')
s <- smc[[1]] + xlim(c(0,400)) + scale_color_manual(values = c('#00bfc4'))+labs(caption = NULL)+ theme_bw(base_size = 16) + theme(legend.position = 'none')
vs <- vaxSMC[[1]] + xlim(c(0,400)) + scale_color_manual(values = c('#c77cff'))+labs(caption = NULL)+ theme_bw(base_size = 16) + theme(legend.position = 'none')
library(patchwork)
plots <- n + v + s + vs + plot_annotation(tag_levels = 'A')
# plots
saveRDS(nothing[[6]], file = paste0('R:/Kelly/synergy_orderly/figures/data/none_trajectories', Sys.Date(), '.rds'))
saveRDS(smc[[6]], file = paste0('R:/Kelly/synergy_orderly/figures/data/smc_trajectories', Sys.Date(), '.rds'))
saveRDS(vax[[6]], file = paste0('R:/Kelly/synergy_orderly/figures/data/vax_trajectories', Sys.Date(), '.rds'))
saveRDS(vaxSMC[[6]], file = paste0('R:/Kelly/synergy_orderly/figures/data/vaxSMC_trajectories', Sys.Date(), '.rds'))
ggsave(filename = 'R:/Kelly/synergy_orderly/figures/trajectories.pdf', plots, height = 7, width = 11, dpi = 500)















# Deterministic version ----
# n_particles = 1
# nothing_det <- run_model(n_particles = n_particles,
#                          n_threads = 1L,
#                          PEV_on = 0,
#                          SMC_on = 0,
#                          tt = 30,
#                          det_mode = TRUE)%>%
#   format_data() %>%
#   make_plots()
# nothing_det[[1]]
# nothing_det[[2]]
# nothing_det[[3]]
# nothing_det[[4]]
# nothing_detplt <- plot_grid(nothing_det[[1]] + labs(title = "No interventions"), nothing_det[[2]], nothing_det[[3]], nothing_det[[4]])
# 
# vax_det <- run_model(n_particles = n_particles,
#                          n_threads = 1L,
#                          PEV_on = 1,
#                          SMC_on = 0,
#                          tt= 30,
#                          det_mode = TRUE,
#                      SMC_kill_rate = 2.5)%>%
#   format_data() %>%
#   make_plots()
# vax_det[[1]]
# vax_det[[2]]
# vax_det[[3]]
# vax_det[[4]]
# vax_detplt <- plot_grid(vax_det[[1]] + labs(title = "Vaccination only"), vax_det[[2]], vax_det[[3]], vax_det[[4]])
# 
# smc_det <- run_model(n_particles = n_particles,
#                          n_threads = 1L,
#                          PEV_on = 0,
#                          SMC_on = 1,
#                          tt= 30,
#                          det_mode = TRUE,
#                      SMC_kill_rate = 2.5)%>%
#   format_data() %>%
#   make_plots()
# smc_det[[1]]
# smc_det[[2]]
# smc_det[[3]]
# smc_det[[4]]
# smc_detplt <- plot_grid(smc_det[[1]] + labs(title = "SMC only"), smc_det[[2]], smc_det[[3]], smc_det[[4]])
# 
# 
# vaxsmc_det <- run_model(n_particles = n_particles,
#                          n_threads = 1L,
#                          PEV_on = 1,
#                          SMC_on = 1,
#                          tt= 30,
#                          det_mode = TRUE,
#                         SMC_kill_rate = 2.5)%>%
#   format_data() %>%
#   make_plots()
# vaxsmc_det[[1]]
# vaxsmc_det[[2]]
# vaxsmc_det[[3]]
# vaxsmc_det[[4]]
# vaxsmc_detplt <- plot_grid(vaxsmc_det[[1]] + labs(title = "Vaccination + SMC"), vaxsmc_det[[2]], vaxsmc_det[[3]], vaxsmc_det[[4]])
# 
# ggsave(filename = "src/plots/nointerventions_det.png", plot = nothing_detplt, width = 11, height = 7)
# ggsave(filename = "src/plots/vaccine_det.png", plot = vax_detplt, width = 11, height = 7)
# ggsave(filename = "src/plots/smc_det.png", plot = smc_detplt, width = 11, height = 7)
# ggsave(filename = "src/plots/vaxsmc_det.png", plot = vaxsmc_detplt, width = 11, height = 7)
# 
# 
# 
# 
# ttoinf <- lapply(list(dfnothing, dfsmc, dfvax, dfvaxsmc), get_ttoinf)
# df_ttoinf <- bind_rows(ttoinf)
# 
# ggplot(df_ttoinf ) + 
#   geom_boxplot(aes(x = scen, y = time*2, color = scen, fill = scen), alpha = 0.3) + 
#   geom_jitter(aes(x = scen, y = time*2, color = scen)) + 
#   theme_bw(base_size = 14) +
#   theme(legend.position = 'none') +
#   scale_y_continuous(breaks = seq(0,max(df_ttoinf$time, na.rm = TRUE)*2, 2)) +
#   labs(x = 'Scenario',
#        y = 'Time to detectable infection (days)') 
# 
# 
# # Plot the SMC efficacy curveo ver time 
# SMC_decay=0.05
# dt = 1:150
# 
# SMC_kill_rate <- 2.3 * exp(-SMC_decay * dt) # time-varying SMC kill rate (exponential decay)
# p_killSMC <- 1 - exp(-SMC_kill_rate * dt) 
# 
# plot(SMC_kill_rate)
# plot(p_killSMC)

