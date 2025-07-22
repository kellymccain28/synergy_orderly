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
                   files = "surv_analysis_model.rds")

tidy_results_trial <- readRDS("surv_analysis_trial.rds")
tidy_results_model <- readRDS('surv_analysis_model.rds')

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
  scale_y_continuous(breaks = seq(-0.3, 1, 0.1),
                     limits = c(-0.3,1),
                     labels = scales::percent) + 
  theme_bw(base_size = 16)

ggsave(filename = 'eff_comparison_model_trial.png', plot = compare)
