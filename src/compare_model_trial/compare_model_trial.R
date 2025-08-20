library(ggplot2)
library(cyphr)
library(tidyverse)
library(orderly2)

orderly_strict_mode()

orderly_dependency(name = 'trial_results',
                   "latest()",
                   files = c('surv_analysis_trial.rds',
                             'monthly_incidence_trial.rds'))
orderly_dependency(name = 'process_model_output',
                   "latest()",
                   files = c('outputs/'))
                   # files = c("outputs/surv_analysis_model.rds",
                   #           "outputs/expected_efficacies.rds",
                   #           "outputs/"))

all_model_outputs <- readRDS('outputs/model_outputs_formatted.rds')

parameters_ll <- readRDS('outputs/parameters_ll.rds')

# Get Sim_id with highest likelihood value 
ml <- parameters_ll[parameters_ll$ll == max(parameters_ll$ll),] 


# Get all simulation directories
output_dir <- 'outputs'
sim_dirs <- list.dirs(output_dir, recursive = FALSE, full.names = TRUE)
sim_dirs <- sim_dirs[grepl("^sim_", basename(sim_dirs))]
sim_dir <- sim_dirs[sim_dirs==paste0('outputs/',ml$ll)]


tidy_results_trial <- readRDS("surv_analysis_trial.rds")
tidy_results_model <- readRDS(paste0(sim_dir, '/surv_analysis_model.rds'))
ve_comparison <- readRDS(paste0(sim_dir, '/expected_efficacies.rds'))

all_results <- rbind(tidy_results_trial %>% mutate(type = 'trial'),
                     tidy_results_model %>% mutate(type = 'model')) %>%
  filter(term %in% c("Both vs. RTSS",'RTSS vs. SMC','Both vs. SMC'))

compare <- ggplot(all_results)+
  geom_point(aes(x = term, y = VE, group = type, color = type), #color = 'darkgreen', 
             position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper, group = type, color = type), 
                width = 0.2, #color = 'darkgreen', 
                position = position_dodge(width = 0.3),
                linewidth = 1) +
  scale_color_manual(values = c('darkgreen','orange')) +
  geom_hline(aes(yintercept = 0)) +
  labs(y = "Efficacy",
       x = NULL,
       color = NULL,
       caption = 'Green is model, red is trial') +
  scale_y_continuous(breaks = seq(-1, 1, 0.1),
                     # limits = c(-0.5,1),
                     labels = scales::percent) + 
  theme_bw(base_size = 12) + 
  facet_wrap(~ year)
# compare
ggsave(filename = 'eff_comparison_model_trial.png', plot = compare, 
       height = 8, width = 8, bg = 'white')


# Compare expected and observed efficacy in model with observed efficacy in trial 
ve_comp_long <- ve_comparison %>%
  select(-c(2:4)) %>%
  pivot_longer(2:4,
               names_to = 'difference_type',
               values_to = 'value') %>%
  mutate(term = 'Both vs. none')

compareobs <- ggplot(all_results)+
  geom_point(aes(x = term, y = VE, group = type, color = type), #color = 'darkgreen', 
             position = position_dodge(width = 0.3), alpha = 0.7) +
  geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper, group = type, color = type), 
                width = 0.2, #color = 'darkgreen', 
                position = position_dodge(width = 0.3),
                linewidth = 1) +
  geom_point(data = ve_comp_long %>% filter(difference_type != 'difference'), 
             aes(x = term, y = value, color = difference_type)) +
  # scale_color_manual(values = c('darkgreen','orange')) +
  geom_hline(aes(yintercept = 0)) +
  labs(y = "Efficacy",
       x = NULL,
       color = NULL,
       caption = 'Green is model, red is trial') +
  scale_y_continuous(breaks = seq(-1, 1, 0.1),
                     # limits = c(-0.5,1),
                     labels = scales::percent) + 
  theme_bw(base_size = 12) + 
  facet_wrap(~ year)
ggsave(filename = 'eff_comparison_model_trial_obs_exp.png', plot = compareobs, 
       height = 8, width = 8, bg = 'white')



# Get the monthly incidence to compare 
modelinci <- readRDS(paste0(sim_dir,'/monthly_incidence_model.rds')) %>%
  # filter(arm !='none') %>%
  mutate(arm = case_when(
    arm == 'vax' ~ 'rtss',
    arm == 'vaxsmc' ~ 'both',
    TRUE ~ arm
  ),
  type = 'model') 
trialinci <- readRDS('monthly_incidence_trial.rds') %>%
  mutate(type = 'trial')

monthlyincidence <- rbind(modelinci, trialinci)

monthlyincidenceplot <- ggplot(data = monthlyincidence) +
  geom_point(aes(x = date, y = incidence_per_1000pm, color = type)) +
  geom_errorbar(aes(x = date, ymin = lower_per_1000, ymax = upper_per_1000, color = type),
                alpha = 0.9, width = 10, linewidth = 0.8) +
  geom_line(aes(x = date, y = incidence_per_1000pm, color = type), 
            linewidth = 0.8, alpha = 0.7) +
  facet_wrap(~arm, nrow = 4,
             scales = 'free_y') +
  scale_y_continuous(breaks = seq(0,150,25)) +
  scale_x_date(breaks = '3 months',
               labels = scales::label_date_short()) + 
  scale_color_manual(values = c('model' = 'lightgreen',
                               'trial' = 'navy')) +
  labs(
    title = "Monthly malaria incidence per 1000 person-months",
    x = "Month",
    y = "Incidence (per 1000 person-months)",
    color = "Study Arm",
    fill = "Study Arm",
    caption = str_glue("Model parameters: \n,
                        SMC max: {ml$max_SMC_kill_rate},\n
                        SMC decay: {ml$SMC_decay},\n
                        season start day: {ml$season_start_day}")
  ) +
  theme_minimal(base_size = 16)

ggsave("trial_monthlyincidence.png", plot = monthlyincidenceplot, bg = 'white')

