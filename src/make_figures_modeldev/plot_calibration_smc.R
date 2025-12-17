# Script to process output of run_smc_test.R
library(ggplot2)
library(tidyverse)

path = "R:/Kelly/synergy_orderly"

# Observed efficacy extracted from figure in SI of Hayley's paper 10.1016/S2214-109X(22)00416-8
observed_efficacy <- read.csv(paste0(path, '/shared/smc_fits_hayley.csv')) %>%
  mutate(weeks_since_smc = ceiling(day_since_smc / 7)) %>%
  group_by(weeks_since_smc) %>%
  mutate(efficacy_week = mean(efficacy))

# this dataset is from a re-run of 3 sets of the best rounds from 2 previous sets of Latin Hypercube sampling with mean least squares 
# ("R:/Kelly/synergy_orderly/src/fit_smc/outputs/grid_search_outputs2025-11-29.rds") and maybe the one from 12-02 ?
# used from test_best_part22025-12-01.rds and test_best2025-12-01
gridsearch <- readRDS('R:/Kelly/synergy_orderly/src/fit_smc/outputs/grid_search2025-12-03.rds') # this one has the best one -- #67 2.333, 16.667, 0.222
gridsearch <- readRDS('R:/Kelly/synergy_orderly/src/fit_smc/outputs/runs_best_smcpars2025-12-03_1.rds')
gridsearch <- readRDS('R:/Kelly/synergy_orderly/src/fit_smc/outputs/runs_best_smcpars2025-12-03_7days.rds')
gridsearch <- readRDS('R:/Kelly/synergy_orderly/src/fit_smc/outputs/runs_best_smcpars_shorterliverstage2025-12-09.rds') # this is with 6 days (?)
gridsearch <- readRDS('R:/Kelly/synergy_orderly/src/fit_smc/outputs/runs_best_smcpars_shorterliverstage2025-12-10.rds') # this is with 8 days 
gridsearch <- readRDS('R:/Kelly/synergy_orderly/src/fit_smc/outputs/runs_best_smcpars_7dayliverstage2025-12-15.rds')
effcumulweekly <- purrr::map_df(gridsearch, 'efficacy',.id = 'sim_id')
pars <- purrr::map_df(gridsearch, 'params')
mlses <- unlist(purrr::map(gridsearch, 'mls',.id = 'sim_id'))
top5 <- order(mlses)[1:2]
mlsesmin <- which(mlses == min(mlses, na.rm = TRUE))
all <- map_df(gridsearch[top5], 'params', .id = 'sim_id')

map_df(gridsearch[mlsesmin], 'params')

# effcumulweekly <- effcumulweekly %>%
#   left_join(pars, by = "sim_id")
# ggplot(effcumulweekly %>% filter(weeks_since_smc < 10, repnum %in% c(1, 7,11,2,16,3,6) &sim_id == 67)) + # for the grid search 12-03
#   geom_point(aes(x = weeks_since_smc, y = efficacy, group = as.factor(sim_id), color = as.factor(repnum)), alpha = 0.4) +
#   geom_line(aes(x = weeks_since_smc, y = efficacy, group = as.factor(sim_id), color = as.factor(repnum)), alpha = 0.3) +
#   geom_hline(aes(yintercept = 0), color = 'darkred', linetype = 2) +
#   geom_line(data = observed_efficacy, aes(x = weeks_since_smc, y = efficacy_week)) +
#   geom_point(data = observed_efficacy, aes(x = weeks_since_smc, y = efficacy_week)) +
#   scale_x_continuous(breaks = seq(0,9,1)) +
#   theme_minimal() +
#   labs(x = 'Weeks since SMC') +
#   facet_wrap(~ repnum)

# if no repnum because just a repetition:
ggplot(effcumulweekly %>% filter(weeks_since_smc < 10) ) + # for the grid search 12-03
  # geom_point(aes(x = weeks_since_smc, y = efficacy, group = as.factor(sim_id)), color = '#709176', alpha = 0.3) +
  geom_line(aes(x = weeks_since_smc, y = efficacy, group = as.factor(sim_id)), color = '#709176', alpha = 0.2) +
  geom_hline(aes(yintercept = 0), color = 'darkred', linetype = 2, linewidth = 1) +
  geom_line(data = observed_efficacy, aes(x = weeks_since_smc, y = efficacy_week), linewidth = 1) +
  geom_point(data = observed_efficacy, aes(x = weeks_since_smc, y = efficacy_week), size = 2) +
  scale_x_continuous(breaks = seq(0,10,1)) +
  scale_y_continuous(labels = scales::label_percent()) +
  theme_minimal(base_size = 14) +
  labs(x = 'Weeks since SMC',
       y = 'Efficacy') + 
  theme(legend.position = 'none')
ggsave(paste0(path, '/figures/smc_fit.pdf'), plot = last_plot(), height = 5, width = 7) # this is 64 runs

matched_complete <- left_join(effcumulweekly %>%
                                rename(predicted_efficacy = efficacy), 
                              observed_efficacy %>% 
                                distinct(weeks_since_smc, efficacy_week) %>% 
                                rename(observed_efficacy = efficacy_week), by = 'weeks_since_smc') %>%
  filter(weeks_since_smc < 10)
mls <- lapply(unique(effcumulweekly$sim_id), function(x){
  d <- matched_complete %>% filter(sim_id == x)
  mean((d$observed_efficacy - d$predicted_efficacy)^2)
})
min(unlist(mls))
minini <- which(mls == min(unlist(mls)))
map_df(gridsearch[minini], 'params')
