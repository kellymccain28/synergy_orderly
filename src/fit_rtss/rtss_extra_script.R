# Script for visualizing/prepping fitting of RTSS parameters 
# simulation script is fit_rtss.R
library(tidyverse)

path = "R:/Kelly/synergy_orderly"
source(file.path(path, 'shared/rtss.R'))

set.seed(123)

# Prepare the "observed" data to which we will compare the cohort sim model output
# Efficacy against infection from White 2015
nsim <- 500
ts <- seq(1,365*3)
phases <- ifelse(ts < 366, 1,
                 ifelse(ts < 729, 2, 3))
csp <- list()
for(i in 1:nsim){
  csp[[i]] <- antibody_titre(t = ts,
                     phase = phases,
                     peak1 = c(621,0.35) ,
                     peak2 = c(277, 0.35),
                     peak3 = c(277,0.35),
                     duration1 = c(45,16),
                     duration2 = c(591,245),
                     rho1 = c(2.37832, 1.00813),
                     rho2 = c(1.034, 1.027),
                     rho3 = c(1.034, 1.027))
}
csp_df <- data.frame(do.call(rbind, csp))
colnames(csp_df) <- ts
csp_df$sim <- 1:nsim

csp_long <- pivot_longer(csp_df, cols = -sim,
                         names_to = "time", values_to = "titre")
csp_long$time <- as.numeric(csp_long$time)

# Plot
# ggplot(csp_long, aes(x = time, y = titre, group = sim)) +
#   geom_line(alpha = 0.3) +
#   theme_minimal() + scale_y_log10() +
#   labs(x = "Time", y = "Antibody Titre")

ve <- lapply(csp, vaccine_eff) # lapply(csp, p_spz_surv, vmin = 0.05)# or should i actually use the vaccine_eff 
ve_df <- data.frame(do.call(rbind, ve))
colnames(ve_df) <- ts
ve_df$sim <- 1:nsim

ve_long <- pivot_longer(ve_df, cols = -sim,
                         names_to = "time", values_to = "ve")
ve_long$time <- as.numeric(ve_long$time)
ve_long <- ve_long %>% ungroup() %>%
  mutate(weeks_since_rtss = ceiling(time/7), # weeks since RTSS because follow-up began 3 weeks after the 3rd dose 
         months_since_rtss = ceiling(time /(365/30))) %>%
  # filter(weeks_since_rtss < 52) %>%
  group_by(weeks_since_rtss, sim) %>%
  summarize(ve_inf = mean(ve))

# Get summary of all efficacy values over time 
ve_long <- ve_long %>% 
  group_by(weeks_since_rtss) %>%
  mutate(observed_efficacy = median(ve_inf),
            efficacy_lower = quantile(ve_inf, 0.025),
            efficacy_higher = quantile(ve_inf, 0.975))

# Plot
ggplot(ve_long) +
  geom_line(aes(x = weeks_since_rtss, y = ve_inf, group = sim), color = '#96BE8C', alpha = 0.2) +
  geom_line(aes(x = weeks_since_rtss, y = observed_efficacy), color = 'darkgreen', linewidth = 1.5) +
  # geom_errorbar(aes(x = weeks_since_rtss, ymin = efficacy_lower, ymax = efficacy_higher)) +
  theme_minimal(base_size = 14) +
  scale_y_continuous(labels = scales::label_percent(),
                     breaks = seq(0,1,0.1), limits = c(0,1),
                     expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0, 156, 10),
                     expand = c(0,0),
                     limits = c(0,156)) +
  geom_vline(xintercept = c(52, 104), linetype = 2) +
  labs(x = "Weeks since 3rd primary dose of RTSS", y = "Efficacy against infection")
ggsave(paste0(path, '/figures/rtss_observed.pdf'), plot = last_plot()) # this is 64 runs

# saveRDS(ve_long %>% select(weeks_since_rtss, observed_efficacy) %>% distinct(), file = paste0(path, '/src/fit_rtss/observed_rtss_efficacy.rds'))
dd <- readRDS(paste0(path, '/src/fit_rtss/observed_rtss_efficacy.rds'))


# Look at the outputs of optimization 
# optimization_results <- readRDS("R:/Kelly/synergy_orderly/src/fit_rtss/outputs/optimization_results_2025-12-03.rds")
# # optimization_results$`1`$eval_history[[20]]$params_tibble %>% select(alpha_ab, beta_ab)
# 
# best_refined <- optimization_results[[which.max(
#   sapply(optimization_results, function(x) x$mls)
# )]]
# 
# best_refined$mls
# best_refined$starting_point_id
# best_refined$initial_params #9.0928452 13.0471435  0.4141795
# best_refined$final_params #9.0739477 15.0000000  0.4586005 # when allowing lambda to go <15, 9.1078524 13.0757521  0.4277123
# best_refined$convergence#0
# best_refined$n_evaluations #231
# eval_history_combined <- optimization_results %>%
#   # Add an index for each optimization run
#   imap_dfr(~ {
#     .x$eval_history %>%
#       map_dfr(function(eval) {
#         params <- eval$params_tibble 
#         params$negll <- eval$mls
#         params$optimization_run <- .y  # Add run identifier
#         params
#       })
#   })


# Look at grid search results 
# "Observed" efficacy from White model 
observed_efficacy_rtss <- readRDS(paste0(path, '/src/fit_rtss/observed_rtss_efficacy.rds'))
rtsseff <- readRDS("R:/Kelly/synergy_orderly/src/fit_rtss/outputs/efficacy_rtss_2025-12-04_1.rds")
parameters <- readRDS("R:/Kelly/synergy_orderly/src/fit_rtss/outputs/parameters_2025-12-04_1.rds")
# rtsscumul <- rtsscumul %>% left_join(parameters)
rtsseff <- rtsseff %>% left_join(parameters)

matched <- observed_efficacy_rtss %>%
  left_join(rtsseff %>% select(weeks_since_rtss, efficacy, sim_id) %>%
              rename(predicted_efficacy = efficacy), 
            by = 'weeks_since_rtss')

# Remove NAs before calculating likelihood
matched_complete <- matched %>%
  filter(!is.na(observed_efficacy), !is.na(predicted_efficacy), !is.infinite(predicted_efficacy))

# Calculate mean least squares for each sim_id 
allmls <- lapply(unique(matched_complete$sim_id),
       function(x){
         m <- matched_complete %>% filter(sim_id == x)
         mean((m$observed_efficacy - m$predicted_efficacy)^2, na.rm = TRUE)
       })
allmls <- unlist(allmls)
which(allmls == min(allmls))
parameters$mls <- allmls


# rtssinfs <- readRDS("R:/Kelly/synergy_orderly/src/fit_rtss/outputs/infectionrecords_rtss_2025-12-02.rds")
best <- rtsseff #%>% filter(alpha_ab > 1.22 & alpha_ab < 1.5)#filter(sim_id %in% parameters[parameters$mls < 0.02,]$sim_id)#
ggplot(best) + 
  geom_line(aes(x = weeks_since_rtss -3, y = efficacy, group = sim_id, color = beta_ab), alpha = 0.3) + 
  # geom_smooth(aes(x = weeks_since_rtss -3, y = efficacy, group = sim_id, color = sim_id)) + 
  geom_line(data = observed_efficacy_rtss, aes(x = weeks_since_rtss, y = observed_efficacy), color = 'black', linewidth= 1) + 
  xlim(c(0,10)) + ylim(c(-0,1)) + 
  theme_minimal() #+
  # theme(legend.position = 'none') 

dr <-p_spz_surv(ab_summary$median_ab, beta_ab = 3.538125, alpha_ab= 1.3804900)
plot(1-dr)
ve <- vaccine_eff(ab_summary$median_ab)
plot(ve)


# Test different parameter values 
library(lhs)
param_ranges <- list(
  alpha_ab = c(0.2, 2),
  beta_ab = c(2,10)
)
A <- randomLHS(100, 2)
# Scale to parameter ranges
params_df <- data.frame(
  alpha_ab = qunif(A[,1], param_ranges$alpha_ab[1], param_ranges$alpha_ab[2]),
  beta_ab = qunif(A[,2], param_ranges$beta_ab[1], param_ranges$beta_ab[2])
)

ts <- seq(0,365*5)
phases <- ifelse(ts < 365, 1,
                 ifelse(ts < 730, 2, 3))
ablist <- map(1:1000, function(i){
  ab_values <- antibody_titre(
    t = ts,
    phase = phases,
    peak1 = c(621, 0.35),
    peak2 = c(277, 0.35),
    peak3 = c(277, 0.35),
    duration1 = c(45, 16),
    duration2 = c(591, 245),
    rho1 = c(2.37832, 1.00813),
    rho2 = c(1.034, 1.027),
    rho3 = c(1.034, 1.027),
    t_boost1 = 365,
    t_boost2 = 730
  )
  
  # Convert to data frame
  tibble(
    time = ts,  # or whatever your time variable is
    ab = ab_values
  )})
# Combine and summarize
ab_summary <- bind_rows(ablist, .id = "rep") %>%
  group_by(time) %>%  # Adjust based on your data structure
  summarise(
    median_ab = median(ab),
    lower_95 = quantile(ab, 0.025),
    upper_95 = quantile(ab, 0.975)
  )

calculate_dr_curves_nested <- function(params_df, ab_df) {
  params_df %>%
    mutate(param_set = row_number()) %>%
    rowwise() %>%
    mutate(
      dr_curve = list(
        ab_df %>%
          mutate(
            DR = 1 - 1 / (1 + (median_ab/ beta_ab)^alpha_ab)
          )
      )
    ) %>%
    ungroup() %>%
    select(param_set, alpha_ab, beta_ab, dr_curve) %>%
    unnest(dr_curve)
}

# Calculate DR curves
dr_curves <- calculate_dr_curves_nested(params_df, ab_summary)

# Plot some examples
dr_curves %>%
  filter(alpha_ab < 1, beta_ab > 6) %>%  
  ggplot(aes(x = time, y = DR, group = param_set, color = alpha_ab)) +
  geom_line() +
  labs(
    title = "DR Curves for Different Parameter Sets",
    x = "Time",
    y = "DR"#,
    # color = "Parameter Set"
  ) +
  theme_minimal() + xlim(c(0,364))
