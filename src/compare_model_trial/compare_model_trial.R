library(ggplot2)
library(cyphr)
library(tidyverse)
library(orderly2)

orderly_strict_mode()

orderly_dependency(name = 'trial_results',
                   "latest()",
                   files = 'surv_analysis_trial.rds')
orderly_dependency(name = 'process_model_output',
                   "latest()",
                   files = c("outputs/surv_analysis_model.rds",
                             "outputs/expected_efficacies.rds"))

tidy_results_trial <- readRDS("surv_analysis_trial.rds")
tidy_results_model <- readRDS('outputs/surv_analysis_model.rds')
ve_comparison <- readRDS('outputs/expected_efficacies.rds')

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

ggsave(filename = 'eff_comparison_model_trial.png', plot = compare, height = 8, width = 8)


# Compare expected and observed efficacy in model with observed efficacy in trial 
ve_comp_long <- ve_comparison %>%
  select(-c(2:4)) %>%
  pivot_longer(2:4,
               names_to = 'difference_type',
               values_to = 'value') %>%
  mutate(term = 'Both vs. none')

compareobs <- ggplot(all_results)+
  geom_point(aes(x = term, y = VE, group = type, color = type), #color = 'darkgreen', 
             position = position_dodge(width = 0.3)) +
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
ggsave(filename = 'eff_comparison_model_trial_obs_exp.png', plot = compareobs, height = 8, width = 8)
