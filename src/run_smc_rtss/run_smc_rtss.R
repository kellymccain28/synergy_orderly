library(tidyverse)
library(ggplot2)
library(odin2)
library(dust2)
library(cowplot)
library(glue)
library(orderly2)

orderly_strict_mode()
orderly_parameters(n_particles = NULL,
                   n_threads = NULL,
                   t_inf = NULL,
                   ts = NULL,
                   tstep = NULL,
                   max_SMC_kill_rate = NULL,
                   lambda = NULL, 
                   kappa = NULL, 
                   # SMC_decay = NULL,
                   season_start_day = NULL,
                   # season_length = NULL,
                   # smc_interval = NULL,
                   inf_start = NULL)

# PEV_on <- 0
# SMC_on <- 0
# n_particles = 1L
# n_threads = 1L
# t_inf = 10#sample(1:(365*3), size = 1) # time of infection relative to vaccination 
# ts = 365 # timesteps (each one is 2 days if ts = 1 below)
VB = 1e6
# det_mode = FALSE
# tstep <- 1
tt <- seq(1, ts, by = tstep)#
threshold <- 100
# max_SMC_kill_rate = 4
# SMC_decay <- 0.05
# SMC_timing <- -10
# SMC parameters
# season_start_day <- 0
# season_length <- 120
# smc_interval <- 30

orderly_shared_resource("smc.R",
                        "rtss.R",
                        "smc_rtss.R",
                        "helper_functions.R")

gen_bs <- odin2::odin("smc_rtss.R")
source("rtss.R")
source("helper_functions.R")

# set.seed(1234)

# Stochastic runs ----
# Without vaccination nor SMC
# n_particles = 500L
# inf_start = 0-- this is when the infection shouldl start relative to the start of follow-up time (more relevant for cohort), can change timing of SMC
nothing <- run_model(n_particles = n_particles,
                     n_threads = 4L,
                     PEV_on = 0,
                     SMC_on = 0,
                     tt= tt,
                     SMC_time = seq(0,100,1),
                     SMC_kill_vec = rep(0,101),
                     t_inf = t_inf,
                     infection_start_day = inf_start,
                     VB = VB,
                     det_mode = FALSE) %>%
  format_data(tt= tt,
              infection_start_day = inf_start,
              n_particles = n_particles) %>%
  make_plots() 
prbc <- nothing[[1]] #+ xlim(c(0,100))
innate <- nothing[[2]]
genad <- nothing[[3]]
varspec <- nothing[[4]]
growthr <- nothing[[5]]
# nothing[[6]]
mplot <- nothing[[7]]
# nothing[[9]]
plot_grid(prbc+labs(caption = NULL, x = 'Days'),
          growthr+labs(caption = NULL, x = 'Days'),
          innate+labs(caption = NULL, x = 'Days'),
          genad+labs(caption = NULL, x = 'Days'),
          varspec+labs(caption = NULL, x = 'Days'),
          mplot+labs(caption = NULL, x = 'Days')
) # in caption, days since start of blood-stage 
ggsave('immunity_plot.pdf')
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

# With vaccination 
vax <- run_model(n_particles = n_particles,
                 n_threads = 4L,
                 PEV_on = 1,
                 SMC_on = 0,
                 t_inf = t_inf,
                 infection_start_day = inf_start,
                 SMC_time = seq(0,100,1),
                 SMC_kill_vec = rep(0,101),
                 tt= tt,
                 VB = VB,
                 det_mode = FALSE) %>%
  format_data(tt= tt,
              infection_start_day = inf_start,
              n_particles = n_particles) %>%
  make_plots()
vax[[1]] + xlim(c(0,100))
vax[[2]]
vax[[3]]
vax[[4]]
vax[[7]]
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
dfvax <- vax[[5]] %>%
  mutate(scen = 'vax')
# % of runs with without at day 60 
# dfvax %>% filter(time == max(time), parasites < 10) %>% count() / n_particles
# dfvax %>% filter(time == 15, parasites < 10) %>% count() / n_particles

# Calculate SMC kill vector 
smc_dose_days <- c(seq(season_start_day, season_start_day + 120 - 1, 30),
                   seq(season_start_day + 365, season_start_day + 365 + 120 - 1, 30),
                   seq(season_start_day + 365*2, season_start_day  + 365*2 + 120 - 1, 30))

time_since_smc <- lapply(list(smc_dose_days), 
                         calc_time_since_dose, 
                         days = 0:(ts*2 -1))
smckillvec <- lapply(time_since_smc,
                                 function(.){
                                   kill <- max_SMC_kill_rate * exp(-(./ lambda)^kappa)  # calculate kill rate with hill function
                                   # kill <- ifelse(is.na(kill), 0, kill)                          
                                   kill[is.na(kill)] <- 0                                       # change NAs (when there is no SMC) to 0 
                                   kill <- kill[seq_along(kill) %% 2 == 1]                      # keep only odd indices since timesteps in the within-host model are 2 days
                                 })

# Without vaccination but with SMC
smc <- run_model(n_particles = n_particles,
                 n_threads = 4L,
                 PEV_on = 0,
                 SMC_on = 1,
                 tt= tt,
                 det_mode = FALSE,
                 t_inf = t_inf,
                 VB = VB,
                 SMC_time = seq(0,length(smckillvec[[1]])-1,1),
                 SMC_kill_vec = smckillvec#,
                 # max_SMC_kill_rate = max_SMC_kill_rate,
                 # SMC parameters
                 # infection_start_day = inf_start,
                 # season_start_day = season_start_day,
                 # season_length = season_length,
                 # smc_interval = smc_interval# smc_timing = 10
                 ) %>%
  format_data(tt= tt,
              infection_start_day = inf_start,
              n_particles = n_particles) %>%
  make_plots()
smc[[1]] #+ xlim(c(0,100))
smc[[2]] #+ xlim(c(0,100))
smc[[3]]
smc[[4]]
# smc[[6]]
smc[[7]]# + labs(x = 'Years', caption = NULL) #+ xlim(c(0,1))
smc[[8]]
smc[[9]]
smc[[10]]
smcplt <- plot_grid(smc[[1]] +labs(caption = NULL, x = 'Days since start of blood stage'), 
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
# % of runs with infection at day 60 
# dfsmc %>% filter(time == max(time), parasites < 10) %>% count() / n_particles
# dfsmc %>% filter(time == 15, parasites < 10) %>% count() / n_particles

# Both vaccination and SMC
vaxSMC <- run_model(n_particles = n_particles,
                    n_threads = 4L,
                    PEV_on = 1,
                    SMC_on = 1,
                    t_inf = t_inf,
                    tt = tt,
                    det_mode = FALSE,
                    VB = VB,
                    SMC_time = seq(0,100,1),
                    SMC_kill_vec = rep(0,101)
                    # max_SMC_kill_rate = max_SMC_kill_rate,
                    # SMC parameters
                    # infection_start_day = inf_start,
                    # season_start_day = season_start_day,
                    # season_length = season_length,
                    # smc_interval = smc_interval# smc_timing = 10
                    ) %>%
  format_data(tt= tt,
              infection_start_day = inf_start,
              n_particles = n_particles) %>%
  make_plots() 
# vaxSMC[[1]]# + xlim(c(0,100))
# vaxSMC[[2]]
# vaxSMC[[3]]
# vaxSMC[[4]]
# vaxSMC[[7]]
# vaxSMC[[9]]
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

saveRDS(dfnothing, file = str_glue("df_nothing_", n_particles,".rds"))
saveRDS(dfvax, file = str_glue("df_vax_", n_particles,".rds"))
saveRDS(dfsmc, file = str_glue("df_smc_", n_particles,".rds"))
saveRDS(dfvaxsmc, file = str_glue("df_vaxsmc_", n_particles,".rds"))

# poster
n <- nothing[[1]] + xlim(c(0,100)) + scale_color_manual(values = c("#F8766D")) +labs(caption = NULL) + theme_bw(base_size = 16) + theme(legend.position = 'none')
v <- vax[[1]] + xlim(c(0,100)) + scale_color_manual(values = c('#7CAE00'))+labs(caption = NULL)+ theme_bw(base_size = 16) + theme(legend.position = 'none')
s <- smc[[1]] + xlim(c(0,100)) + scale_color_manual(values = c('#00bfc4'))+labs(caption = NULL)+ theme_bw(base_size = 16) + theme(legend.position = 'none')
vs <- vaxSMC[[1]] + xlim(c(0,100)) + scale_color_manual(values = c('#c77cff'))+labs(caption = NULL)+ theme_bw(base_size = 16) + theme(legend.position = 'none')
library(patchwork)
plots <- n + v + s + vs
# plots
ggsave(filename = 'trajectories.png', plots, height = 5.5, width = 11)















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

