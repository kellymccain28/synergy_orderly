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
lhs_parameters <- readRDS("R:/Kelly/synergy_orderly/src/fit_smc/outputs/test_fitted_params_smc_2025-11-19.rds")#lhs_parameters_0411.rds
effweekly <- purrr::map_df(lhs_parameters, 'efficacy_weekly')
effdaily <- purrr::map_df(lhs_parameters, 'efficacy_daily')
effcumulweekly <- purrr::map_df(lhs_parameters, 'efficacy_weekly_cumul')
effcumuldaily <- purrr::map_df(lhs_parameters, 'efficacy_daily_cumul')
lhspars <- purrr::map_df(lhs_parameters, 'params')
lhsparasit <- purrr::map_df(lhs_parameters, 'parasitemia', .id = "sim_id")
lhsinfrecords <- purrr::map_df(lhs_parameters, 'infection_records', .id = 'sim_id')

# effcumulweekly <- calc_smc_efficacy_cumul(lhsinfrecords, params_row = lhspars[1,])

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
ggplot() + 
  geom_line(data = observed_efficacy, aes(x = ceiling(day_since_smc/7), y = efficacy), color = 'orchid') +
  geom_line(data = observed_efficacy, aes(x = weeks_since_smc, y = efficacy_week), color = 'orchid4') +
  geom_line(data = observed_efficacy_week, aes(x = weeks_since_smc, y = observed_efficacy), color = 'orange') + 
  theme_bw()
ggplot() + 
  geom_line(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), color = 'orchid') +
  geom_line(data = observed_efficacy_week, aes(x = weeks_since_smc *7+3, y = observed_efficacy), color = 'orange') + 
  geom_point(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), color = 'orchid') +
  geom_point(data = observed_efficacy_week, aes(x = weeks_since_smc *7+3, y = observed_efficacy), color = 'orange') + 
  theme_bw()
ggplot(observed_efficacy) + 
  geom_line(aes(x = ceiling(day_since_smc/7), y = efficacy), color = 'orchid4') +
  geom_point(aes(x = weeks_since_smc, y = eff_week), color = 'orange') + 
  geom_line(aes(x = day_since_smc, y = efficacy), color = 'orchid4') +
  geom_point(aes(x = weeks_since_smc*7, y = eff_week), color = 'darkorange') +
  theme_bw()
ggplot()+
  geom_point(data =observed_efficacy, aes(x = weeks_since_smc, y = eff_week)) + 
  geom_point(data  =observed_efficacy_week, aes(x = weeks_since_smc+1, y = observed_efficacy), color = 'blue')

# # Look at efficacy output when summarizing by cumulative proportion of population infected (not total incidence)
ggplot(effcumuldaily %>% filter(days_since_smc < 70)) +
  geom_point(aes(x = days_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.2) +
  geom_line(aes(x = days_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.2) +
  ylim(c(-0.5, 1)) + #xlim(c(115, 170)) +
  geom_line(data = observed_efficacy, aes(x = day_since_smc, y = efficacy)) +
  geom_point(data = observed_efficacy, aes(x = day_since_smc, y = efficacy)) +
  theme_minimal() +
  theme(legend.position = 'none')

ggplot(effcumulweekly %>% filter(weeks_since_smc < 10)) +
  # geom_point(data = effcumuldaily%>% filter(days_since_smc < 70), aes(x = days_since_smc/7, y = efficacy, group = sim_id), color = 'grey', alpha = 0.1) +
  # geom_line(data = effcumuldaily%>% filter(days_since_smc < 70), aes(x = days_since_smc/7, y = efficacy, group = sim_id), color = 'grey', alpha = 0.1) +
  geom_point(aes(x = weeks_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.4) +
  geom_line(aes(x = weeks_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.3) +
  ylim(c(-0.3, 1)) +
  geom_hline(aes(yintercept = 0), color = 'darkred', linetype = 2) +
  geom_line(data = observed_efficacy, aes(x = weeks_since_smc, y = efficacy_week)) +
  geom_point(data = observed_efficacy, aes(x = weeks_since_smc, y = efficacy_week)) +
  geom_line(data = observed_efficacy, aes(x = day_since_smc/7, y = efficacy)) +
  geom_point(data = observed_efficacy, aes(x = day_since_smc/7, y = efficacy)) +
  scale_x_continuous(breaks = seq(0,9,1)) +
  theme_minimal() +
  labs(x = 'Weeks since SMC') +
  theme(legend.position = 'none')

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

ggplot(plot_data, aes(x = time, y = survival)) +
  geom_line(color = "steelblue", linewidth = 1) +
  facet_wrap(~ label) +
  labs(
    title = "Weibull Survival Functions",
    x = "Time",
    y = "Kill rate (t))"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = 'none'
  )

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
optimization_results <- readRDS("R:/Kelly/synergy_orderly/src/fit_smc/outputs/optimization_results_2210.rds") # 2210 does allow superinfections

best_refined <- optimization_results[[which.max(
  sapply(optimization_results, function(x) x$log_likelihood)
)]]

best_refined$log_likelihood
best_refined$starting_point_id
best_refined$initial_params #9.0928452 13.0471435  0.4141795
best_refined$final_params #9.0739477 15.0000000  0.4586005 # when allowing lambda to go <15, 9.1078524 13.0757521  0.4277123
best_refined$convergence#0
best_refined$n_evaluations #231
pars_to_test <- data.frame(
  max_SMC_kill_rate = rep(best_refined$final_params[1], 10),
  lambda = rep(best_refined$final_params[2],10),
  kappa = rep(best_refined$final_params[3],10),
  sim_id = seq(1, 10),
  lag_p_bite = rep(0,10),
  season_start_day = rep(50,10),
  smc_dose_days = rep(10,10),
  p_bite = I(rep(params_list[[1]]$p_bite, 10))
)
pars_to_test <- split(pars_to_test, seq(nrow(pars_to_test)))
results2 <- lapply(pars_to_test,
                   function(params_row){
                     o <- run_cohort_simulation(params_row, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                                                metadata_df,
                                                base_inputs,
                                                output_dir = 'R:/Kelly/src/fit_smc/simulation_outputs',
                                                allow_superinfections = TRUE,
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

                     return(list(efficacy_weekly = eff,
                                 efficacy_daily = eff_daily,
                                 params = params_row,
                                 parasitemia = o$parasitemia_data,
                                 infection_records_smc = o$infection_records))
                   })
eff_weekly <- purrr::map_df(results2, 'efficacy_weekly')
eff_daily <- purrr::map_df(results2, 'efficacy_daily')
ggplot(eff_daily %>% filter(days_since_smc < 70)) +
  geom_point(aes(x = days_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.8) +
  geom_line(aes(x = days_since_smc, y = efficacy, group = sim_id, color = sim_id), alpha = 0.2) +
  ylim(c(-0.5, 1)) +
  geom_line(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), alpha = 0.2) +
  geom_point(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), alpha = 0.2) +
  theme_minimal() +
  theme(legend.position = 'none')
ggplot(eff_weekly %>% filter(weeks_since_smc < 10)) +
  ylim(c(-0.5, 1)) +
  geom_line(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), alpha = 0.8) +
  geom_point(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), alpha = 0.8) +
  geom_point(aes(x = weeks_since_smc*7, y = efficacy, group = sim_id, color = as.factor(sim_id)), alpha = 0.8) +
  geom_line(aes(x = weeks_since_smc*7, y = efficacy, group = sim_id, color = as.factor(sim_id)), alpha = 0.8) +
  theme_minimal() +
  theme(legend.position = 'none')

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

