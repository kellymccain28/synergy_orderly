# Script to process output of run_smc_test.R
library(ggplot2)
library(tidyverse)

path = "R:/Kelly/synergy_orderly"

# Observed efficacy extracted from figure in SI of Hayley's paper 10.1016/S2214-109X(22)00416-8
observed_efficacy <- read.csv(paste0(path, '/shared/smc_fits_hayley.csv')) %>%
  mutate(weeks_since_smc = ceiling(day_since_smc / 7)) %>%
  group_by(weeks_since_smc) %>%
  mutate(efficacy_week = mean(efficacy))
observed_efficacy_week <- read.csv(paste0(path, '/shared/smc_fits_hayley_week.csv'))
lhs_parameters <- readRDS("R:/Kelly/synergy_orderly/src/fit_smc/outputs/test_best_part22025-12-01.rds")##lhs_parameters_0411.rds
lhs_parameters2 <- readRDS("R:/Kelly/synergy_orderly/src/fit_smc/outputs/grid_search2025-12-01.rds")
# lhs_parameters <- readRDS("R:/Kelly/synergy_orderly/src/fit_smc/outputs/grid_search2025-12-02.rds") 
# effweekly <- purrr::map_df(lhs_parameters, 'efficacy_weekly')
# effdaily <- purrr::map_df(lhs_parameters, 'efficacy_daily')
effcumulweekly <- purrr::map_df(lhs_parameters, 'efficacy', .id = "sim_id")
lhspars <- purrr::map_df(lhs_parameters, 'params')
lhsparasit <- purrr::map_df(lhs_parameters, 'parasitemia', .id = "sim_id")
lhsinfrecords <- purrr::map_df(lhs_parameters, 'infection_records', .id = 'sim_id')

# Look at MLS values 
mlses <- unlist(purrr::map(lhsparameters, 'mls',.id = 'sim_id'))
mlses1 <- purrr::map(lhs_parameters1, 'mls', .id = 'sim_id')
mlses1 <- unlist(mlses1)
mlses2 <- purrr::map(lhs_parameters2, 'mls', .id = 'sim_id')
mlses2 <- unlist(mlses2)
mlses <- c(mlses1, mlses2)
mlsesmin <- which(mlses == min(mlses, na.rm = TRUE))
lowestmls1 <- which(mlses1 == min(mlses1, na.rm = TRUE))
lowestmls2 <- which(mlses2 == min(mlses2, na.rm = TRUE))
lhs_parameters1[[lowestmls1]]$params[1:3] # 2.071429, 20.71429, 0.5 from first grid search
lhs_parameters2[[lowestmls2]]$params[1:3] # 2.071429, 20.71429, 0.5 from first grid search
lhsparameters[[mlsesmin]]$params[1:3]
# get top 5 
top5 <- order(mlses)[1:10]
top51 <- order(mlses1)[1:10] # get top 5 
top52 <- order(mlses2)[1:20]
paramstop5 <- map_df(lhsparameters[top5], "params") %>% select(max_SMC_kill_rate, lambda, kappa) %>%
  mutate(mls = mlses[top5], set = 'new')
paramstop51 <- map_df(lhs_parameters1[top51], "params") %>% select(max_SMC_kill_rate, lambda, kappa) %>%
  mutate(mls = mlses1[top51], set = 1)
paramstop52 <- map_df(lhs_parameters2[top52], "params") %>% select(max_SMC_kill_rate, lambda, kappa) %>%
  mutate(mls = mlses2[top52], set = 2)
paramstop10 <- rbind(paramstop51, paramstop52) %>% arrange(mls)

saveRDS(paramstop10, 'R:/Kelly/synergy_orderly/src/fit_smc/outputs/top_20params_1201.rds')

paramsall <- map_df(lhs_parameters, "params") %>% select(max_SMC_kill_rate, lambda, kappa) %>%
  mutate(mls = mlses)

gridsearch <- readRDS('R:/Kelly/synergy_orderly/src/fit_smc/outputs/grid_search2025-12-03.rds') # this one has the best one -- #67 2.333, 16.667, 0.222
effcumulweekly <- purrr::map_df(gridsearch, 'efficacy')
pars <- purrr::map_df(gridsearch, 'params')
mlses <- unlist(purrr::map(gridsearch, 'mls',.id = 'sim_id'))
top5 <- order(mlses)[1:10]
mlsesmin <- which(mlses == min(mlses, na.rm = TRUE))
all <- map_df(gridsearch, 'params')
# # Look at efficacy output
# ggplot(effdaily %>% filter(days_since_smc < 70 & days_since_smc %% 2 == 0)) +
#   geom_point(aes(x = days_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.2) +
#   geom_line(aes(x = days_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.2) +
#   ylim(c(-0.3, 1)) + #xlim(c(115, 170)) +
#   geom_line(data = observed_efficacy, aes(x = day_since_smc, y = efficacy)) +
#   geom_point(data = observed_efficacy, aes(x = day_since_smc, y = efficacy)) +
#   theme_minimal() +
#   theme(legend.position = 'none')
# 
# ggplot(effweekly %>% filter(weeks_since_smc < 10)) +
#   geom_point(aes(x = weeks_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.1) +
#   geom_line(aes(x = weeks_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.2) +
#   ylim(c(-0.5, 1)) +
#   geom_line(data = observed_efficacy, aes(x = day_since_smc/7, y = efficacy)) +
#   geom_point(data = observed_efficacy, aes(x = day_since_smc/7, y = efficacy)) +
#   scale_x_continuous(breaks = seq(0,9,1)) +
#   theme_minimal() +
#   theme(legend.position = 'none')

# Compare weekly and daily summarized observed efficacy 
# ggplot() + 
#   geom_line(data = observed_efficacy, aes(x = ceiling(day_since_smc/7), y = efficacy), color = 'orchid') +
#   geom_line(data = observed_efficacy, aes(x = weeks_since_smc, y = efficacy_week), color = 'orchid4') +
#   geom_line(data = observed_efficacy_week, aes(x = weeks_since_smc, y = observed_efficacy), color = 'orange') + 
#   theme_bw()
# ggplot() + 
#   geom_line(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), color = 'orchid') +
#   geom_line(data = observed_efficacy_week, aes(x = weeks_since_smc *7+3, y = observed_efficacy), color = 'orange') + 
#   geom_point(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), color = 'orchid') +
#   geom_point(data = observed_efficacy_week, aes(x = weeks_since_smc *7+3, y = observed_efficacy), color = 'orange') + 
#   theme_bw()
# ggplot(observed_efficacy) + 
#   geom_line(aes(x = ceiling(day_since_smc/7), y = efficacy), color = 'orchid4') +
#   geom_point(aes(x = weeks_since_smc, y = eff_week), color = 'orange') + 
#   geom_line(aes(x = day_since_smc, y = efficacy), color = 'orchid4') +
#   geom_point(aes(x = weeks_since_smc*7, y = eff_week), color = 'darkorange') +
#   theme_bw()
# ggplot()+
#   geom_point(data =observed_efficacy, aes(x = weeks_since_smc, y = eff_week)) + 
#   geom_point(data  =observed_efficacy_week, aes(x = weeks_since_smc+1, y = observed_efficacy), color = 'blue')

# # Look at efficacy output when summarizing by cumulative proportion of population infected (not total incidence)
# ggplot(effcumuldaily %>% filter(days_since_smc < 70 & sim_id%in% lhspars$sim_id[11:20])) +
#   geom_point(aes(x = days_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.2) +
#   geom_line(aes(x = days_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.2) +
#   ylim(c(-0.5, 1)) + #xlim(c(115, 170)) +
#   geom_line(data = observed_efficacy, aes(x = day_since_smc, y = efficacy)) +
#   geom_point(data = observed_efficacy, aes(x = day_since_smc, y = efficacy)) +
#   theme_minimal() +
#   theme(legend.position = 'none')
effcumulweekly <- effcumulweekly %>%
  left_join(pars, by = "sim_id")
ggplot(effcumulweekly %>% filter(weeks_since_smc < 10, repnum %in% c(1, 7,11,2,16,3,6) &sim_id == 67)) + # for the grid search 12-03
  # geom_point(data = effcumuldaily%>% filter(days_since_smc < 70), aes(x = days_since_smc/7, y = efficacy, group = sim_id), color = 'grey', alpha = 0.1) +
  # geom_line(data = effcumuldaily%>% filter(days_since_smc < 70), aes(x = days_since_smc/7, y = efficacy, group = sim_id), color = 'grey', alpha = 0.1) +
  geom_point(aes(x = weeks_since_smc, y = efficacy, group = as.factor(sim_id), color = as.factor(repnum)), alpha = 0.4) +
  geom_line(aes(x = weeks_since_smc, y = efficacy, group = as.factor(sim_id), color = as.factor(repnum)), alpha = 0.3) +
  # ylim(c(-0.3, 1)) +
  geom_hline(aes(yintercept = 0), color = 'darkred', linetype = 2) +
  geom_line(data = observed_efficacy, aes(x = weeks_since_smc, y = efficacy_week)) +
  geom_point(data = observed_efficacy, aes(x = weeks_since_smc, y = efficacy_week)) +
  # geom_line(data = observed_efficacy, aes(x = day_since_smc/7, y = efficacy)) +
  # geom_point(data = observed_efficacy, aes(x = day_since_smc/7, y = efficacy)) +
  scale_x_continuous(breaks = seq(0,9,1)) +
  theme_minimal() +
  labs(x = 'Weeks since SMC') +
  facet_wrap(~ repnum)
  # theme(legend.position = 'none')
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

# # Assume observed data has Normal measurement error
calculate_likelihood <- function(observed, predicted, sigma = 0.05) {
  # neg log-likelihood
  -sum(dnorm(observed, mean = predicted, sd = sigma, log = TRUE), na.rm = TRUE)
}

eff <- effweekly
eff <- effdaily
likelihoods <- eff %>%
  ungroup() %>%
  mutate(#day_since_smc = round(weeks_since_smc * 7,0),
    day_since_smc = days_since_smc
  ) %>%
  rename(predicted_efficacy = efficacy) %>%
  left_join(observed_efficacy %>% rename(observed_efficacy = efficacy), by = "day_since_smc") %>%
  group_by(sim_id) %>%
  summarise(
    neg_log_lik = calculate_likelihood(observed_efficacy,
                                       predicted_efficacy,
                                       sigma = 0.1),  # Adjust sigma based on your uncertainty
    .groups = 'drop'
  ) %>%
  mutate(
    likelihood = exp(-neg_log_lik),
    weight = likelihood / sum(likelihood)  # Normalize to get weights
  ) %>%
  arrange(neg_log_lik)

# Get top 10 runs
top_runs <- likelihoods %>% slice_head(n = 10)

# Plot them
effcumuldaily %>%
  filter(sim_id %in% top_runs$sim_id &
           days_since_smc < 10*7 &
           days_since_smc %%2 == 0
         # weeks_since_smc < 10
  ) %>%
  # ggplot(aes(x = weeks_since_smc*7, y = efficacy, group = sim_id, color = as.factor(sim_id))) +
  ggplot(aes(x = days_since_smc, y = efficacy, group = sim_id, color = as.factor(sim_id))) +
  geom_line(alpha = 0.5) +
  geom_point(data = observed_efficacy, aes(x = day_since_smc, y = efficacy),
             color = "black", inherit.aes = FALSE) +
  geom_line(data = observed_efficacy, aes(x = day_since_smc, y = efficacy),
            color = "black", inherit.aes = FALSE) +
  labs(title = "Top 10 runs vs observed SMC efficacy",
       y = "SMC efficacy", x = "Days since SMC") +
theme_minimal() + ylim(c(-0.5, 1)) + theme(legend.position = 'none')
#
# # parameter_set_293_generic_TRUE is closest
# # max_SMC_kill_rate lambda kappa season_start_day sim_id                         lag_p_bite
# # <dbl>             <dbl>  <dbl> <dbl>            <chr>                          <dbl>
# # 6.41              44.2   0.425 150              parameter_set_293_generic_TRUE 0
#
# top 10
# max_SMC_kill_rate lambda kappa season_start_day sim_id                         lag_p_bite
# <dbl>  <dbl> <dbl>            <dbl> <chr>                               <dbl>
#   1             12.7   41.2  8.20                98 parameter_set_114_generic_TRUE          0
# 2              8.80  36.8  9.55                51 parameter_set_133_generic_TRUE          0
# 3              7.99   9.99 0.235               41 parameter_set_248_generic_TRUE          0
# 4              6.41  44.2  0.425              150 parameter_set_293_generic_TRUE          0
# 5             24.5   30.0  3.71                99 parameter_set_420_generic_TRUE          0
# 6             20.0   36.1  7.71                89 parameter_set_545_generic_TRUE          0
# 7             21.0   32.1  4.15                37 parameter_set_802_generic_TRUE          0
# 8             16.0   33.5  8.39                74 parameter_set_867_generic_TRUE          0
# 9             14.4   33.5  9.87                94 parameter_set_934_generic_TRUE          0
# 10             22.7   31.0  9.05               123 parameter_set_993_generic_TRUE          0

eff %>%
  filter(sim_id %in% top_runs$sim_id &
           days_since_smc < 10*5 &
           days_since_smc %% 2== 0
         # weeks_since_smc < 10
  ) %>%
  ggplot() +
  # geom_line(aes(x = weeks_since_smc, y = inci_none, group = sim_id, color = as.factor(sim_id)),
  #           alpha = 0.8, linetype = 2) +
  # geom_line(aes(x = weeks_since_smc, y = inci_smc, group = sim_id, color = as.factor(sim_id)),
  #           alpha = 0.8) +
  geom_line(aes(x = days_since_smc, y = inci_none, group = sim_id, color = as.factor(sim_id)),
            alpha = 0.8, linetype = 2) +
  geom_line(aes(x = days_since_smc, y = inci_smc, group = sim_id, color = as.factor(sim_id)),
            alpha = 0.8) +
  labs(y = "incidence", x = "days since SMC") +
  theme_minimal() #+ ylim(c(0,0.7))

weibull_survival <- function(t, max_SMC_kill_rate, lambda, kappa){
  max_SMC_kill_rate * exp(-(t/lambda)^kappa)
}

t_seq = seq(0, 60, length.out = 60)
pars <- lhspars %>% filter(sim_id %in% top_runs$sim_id)

pars <- data.frame(
  max_SMC_kill_rate = rep(3,10),
  lambda = rep(16, 10),
  kappa= seq(0.3,0.7, length.out = 10),
  sim_id = paste0(seq(1,10), 'id')
)
plot_data <- pars %>%
  mutate(numid = readr::parse_number(sim_id)) %>%
  mutate(label = paste0(numid, "\n", round(max_SMC_kill_rate,2), ', ', round(lambda,2),', ', round(kappa,2))) %>%
  rowwise() %>%
  do({
    data.frame(
      id = .$sim_id,
      label = .$label,
      time = t_seq,
      survival = weibull_survival(t_seq, .$max_SMC_kill_rate, .$lambda, .$kappa)
    )
  }) %>%
  ungroup()

ggplot(plot_data, aes(x = time, y = survival, group = label, color = label)) +
  geom_line(linewidth = 0.8) +#color = "steelblue"
  # facet_wrap(~ label) +
  labs(
    title = "Weibull Survival Functions",
    x = "Time",
    y = "Kill rate (t))"
  ) +
  theme_minimal() #+
  # theme(
  #   strip.text = element_text(size = 11, face = "bold"),
  #   plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
  #   legend.position = 'none'
  # )

# Try to understand why the cohort sims are resulting in DROPs in efficacy all of the sudden  -- SEEMS TO BE RELATED TO THE KILL RATES THAT ARE ALL DROPPING PRECIPITOUSLY
test <- o$parasitemia_data %>%#lhsparasit %>%
  # filter(sim_id == 5) %>%
  filter(rid == 5) %>%
  mutate(detectable = ifelse(any(parasites> threshold), 1, 0))
table(test$detection_day)
table(test$day1_BSinfection)
table(test$rid)
table(test$arm)

ggplot(test) + 
  geom_line(aes(x = time_ext, y = prob_smckill, group = paste0(rid, ', ', day1_BSinfection)))+xlim(c(0,50))
ggplot(test)+
  geom_line(aes(x = time_ext, y = parasites, group = day1_BSinfection, color = as.factor(day1_BSinfection),
                linetype = as.factor(detectable))) + xlim(c(0,50))+ 
  scale_y_log10()
ggplot(test)+
  geom_line(aes(x = time_ext, y = SMC_kill_rateout, group = day1_BSinfection, color = as.factor(day1_BSinfection),
                linetype = as.factor(detectable)))  + xlim(c(0,70))
ggplot() + 
  +  geom_line(aes(x = 1:70, y = ratetest*2))



# Optimization
optimization_results <- readRDS("R:/Kelly/synergy_orderly/src/fit_smc/outputs/optimization_results_20251128.rds") # 2210 does allow superinfections

best_refined <- optimization_results[[which.max(
  sapply(optimization_results, function(x) x$mean_least_squares)
)]]

best_refined$mean_least_squares
best_refined$starting_point_id
best_refined$initial_params #9.0928452 13.0471435  0.4141795
best_refined$final_params #9.0739477 15.0000000  0.4586005 # when allowing lambda to go <15, 9.1078524 13.0757521  0.4277123
best_refined$convergence#0
best_refined$n_evaluations #231
reps <- 3
pars_to_test <- data.frame(
  max_SMC_kill_rate = rep(best_refined$initial_params[1], reps),
  lambda = rep(best_refined$final_params[2],reps),
  kappa = rep(best_refined$final_params[3],reps),
  sim_id = seq(1, reps),
  lag_p_bite = rep(0,reps),
  season_start_day = rep(50,reps),
  smc_dose_days = rep(reps,reps),
  p_bite = I(rep(params_list[[1]]$p_bite, reps))
)
eval_history_combined <- optimization_results %>%
  # Add an index for each optimization run
  imap_dfr(~ {
    .x$eval_history %>%
      map_dfr(function(eval) {
        params <- eval$params_tibble
        params$mls <- eval$mls
        params$
        params$optimization_run <- .y  # Add run identifier
        params
      })
  })
final_params_df <- optimization_results %>%
  map_dfr(~ {
    tibble(
      starting_point_id = .x$starting_point_id,
      max_SMC_kill_rate = .x$final_params[1],
      lambda = .x$final_params[2],
      kappa = .x$final_params[3],
      final_mls = .x$mls,
      convergence = .x$convergence,
      n_evaluations = .x$n_evaluations
    )
  })
pars_to_test <- data.frame(
  max_SMC_kill_rate = c(final_params_df$max_SMC_kill_rate, 2.071429), 
  lambda = c(final_params_df$lambda, 20.71429),
  kappa = c(final_params_df$kappa, 0.5),
  sim_id = c(final_params_df$starting_point_id,'grid'),
  lag_p_bite = rep(0,3),
  season_start_day = rep(0,3),
  smc_dose_days = rep(10,3),
  p_bite = I(rep(params_list[[1]]$p_bite, 3))
)
# pars_to_test <- pars_to_test[1,]
pars_to_test <- split(pars_to_test, seq(nrow(pars_to_test)))
results2 <- lapply(pars_to_test,
                   function(params_row){
                     o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                                                metadata_df,
                                                base_inputs,
                                                output_dir = 'R:/Kelly/src/fit_smc/simulation_outputs',
                                                # allow_superinfections = TRUE,
                                                return_parasitemia = TRUE,
                                                save_outputs = FALSE)
                  
                     eff <- calc_smc_efficacy(o$infection_records,
                                              params_row,
                                              by_week = TRUE)
                     eff_daily <- calc_smc_efficacy(o$infection_records,
                                                    params_row,
                                                    by_week = FALSE)
                     eff$sim_id <- params_row$sim_id
                     eff_daily$sim_id <- params_row$sim_id
                     # Efficacy with cumulative proportion
                     eff_cumul <- calc_smc_efficacy_cumul(o$infection_records,
                                                          params_row,
                                                          by_week = TRUE)
                     eff_daily_cumul <- calc_smc_efficacy_cumul(o$infection_records,
                                                                params_row,
                                                                by_week = FALSE)
                     eff_cumul$sim_id <- params_row$sim_id
                     eff_daily_cumul$sim_id <- params_row$sim_id
                     
                     return(list(efficacy_weekly = eff,
                                 efficacy_daily = eff_daily,
                                 efficacy_weekly_cumul = eff_cumul, 
                                 efficacy_daily_cumul = eff_daily_cumul,
                                 parasitemia = o$parasitemia_data %>% mutate(sim_id = params_row$sim_id),
                                 infection_records = o$infection_records %>% mutate(sim_id = params_row$sim_id),
                                 params = params_row))
                   })
eff_weekly <- purrr::map_df(results2, 'efficacy_weekly')
eff_daily <- purrr::map_df(results2, 'efficacy_daily')
effcumulweekly <- purrr::map_df(results2, 'efficacy_weekly_cumul')
infrecords <- purrr::map_df(results2, 'infection_records')

filt <- infrecords %>% filter(time_ext >=0) 
eff1 <- calc_smc_efficacy_cumul(filt %>% filter(sim_id == 1),#o$infection_records,
                                params_row,
                                by_week = TRUE) %>% mutate(sim_id = 1)
eff2 <- calc_smc_efficacy_cumul(filt %>% filter(sim_id == 2),#o$infection_records,
                                params_row,
                                by_week = TRUE) %>% mutate(sim_id = 2)
effgrid <- calc_smc_efficacy_cumul(filt %>% filter(sim_id == 'grid'),#o$infection_records,
                                params_row,
                                by_week = TRUE) %>% mutate(sim_id = 'grid')
effcumulweekly <- rbind(eff1, eff2, effgrid)

# ggplot(eff_daily %>% filter(days_since_smc < 70)) +
#   geom_point(aes(x = days_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.8) +
#   geom_line(aes(x = days_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.2) +
#   ylim(c(-0.5, 1)) +
#   geom_line(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), alpha = 0.2) +
#   geom_point(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), alpha = 0.2) +
#   theme_minimal() +
#   theme(legend.position = 'none')
# ggplot(eff_weekly %>% filter(weeks_since_smc < 10)) +
#   ylim(c(-0.5, 1)) +
#   geom_line(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), alpha = 0.8) +
#   geom_point(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), alpha = 0.8) +
#   geom_point(aes(x = weeks_since_smc*7, y = efficacy, group = sim_id, color = as.factor(sim_id)), alpha = 0.8) +
#   geom_line(aes(x = weeks_since_smc*7, y = efficacy, group = sim_id, color = as.factor(sim_id)), alpha = 0.8) +
#   theme_minimal() +
#   theme(legend.position = 'none')
ggplot(effcumulweekly %>% filter(weeks_since_smc < 10)) +
  # ylim(c(-0.5, 1)) +
  geom_line(data = observed_efficacy, aes(x = weeks_since_smc, y = efficacy), alpha = 0.8) +
  geom_point(data = observed_efficacy, aes(x = weeks_since_smc, y = efficacy), alpha = 0.8) +
  geom_point(aes(x = weeks_since_smc, y = efficacy, color = as.factor(sim_id)),  alpha = 0.8) +
  geom_line(aes(x = weeks_since_smc, y = efficacy, group = as.factor(sim_id), color = as.factor(sim_id)), alpha = 0.8) +
  theme_minimal() +
  theme(legend.position = 'none')
# observed_efficacy %>%
#   mutate(weeks_since_smc = ceiling(day_since_smc / 7)) %>%
#   group_by(weeks_since_smc) %>%
#   summarize(observed_efficacy = mean(efficacy)) %>%
#   ggplot()+geom_line(aes(x = weeks_since_smc, y = observed_efficacy))
ggplot(eff_weekly %>% filter(weeks_since_smc < 10)) +
    geom_line(aes(x = weeks_since_smc, y = inci_none, group = sim_id, color = as.factor(sim_id)),
              alpha = 0.8, linetype = 2) +
    geom_line(aes(x = weeks_since_smc, y = inci_smc, group = sim_id, color = as.factor(sim_id)),
              alpha = 0.8) +
    labs(y = "incidence", x = "weeks since SMC") +
    theme_minimal() +#+ ylim(c(0,0.7))
    theme(legend.position = 'none')

ggplot(eff_daily %>% filter(days_since_smc < 70)) +
  geom_line(aes(x = days_since_smc, y = inci_none, group = sim_id, color = as.factor(sim_id)),
            alpha = 0.8, linetype = 2) +
  geom_line(aes(x = days_since_smc, y = inci_smc, group = sim_id, color = as.factor(sim_id)),
            alpha = 0.8) +
  labs(y = "incidence", x = "days since SMC") +
  theme_minimal() +#+ ylim(c(0,0.7))
  theme(legend.position = 'none')

