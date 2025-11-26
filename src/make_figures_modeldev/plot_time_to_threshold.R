# Script to plot the time it takes for each infection to reach the threshold, by arm 

# Load packages
library(ggplot2)
library(colorspace)

arm_colors <- c('#E1BC29',
                '#7768AE',
                '#3BB273',
                '#E15554')
lighter_arm_colors <- colorspace::lighten(arm_colors, amount = 0.3)

path <- 'R:/Kelly/synergy_orderly/src/sim_cohort_generic/outputs/'
outputsfolder <- 'outputs_2025-11-26'

# Using the outputs from monthly_incidence_plot.R
formatted <- readRDS(paste0(path, outputsfolder, '/formatted_infrecords.rds')) %>%
  filter(detectable==1 ) %>%
  mutate(arm = factor(arm, levels = c('none','smc','rtss','both'))) 


ggplot(formatted) + 
  geom_boxplot(aes(x = arm, y = t_toreach_threshold, color = arm, fill = arm), linewidth = 0.8) + 
  # geom_jitter(aes(x = arm, y = t_toreach_threshold, color = arm), alpha = 0.2) + 
  scale_y_log10() + 
  scale_color_manual(values =  arm_colors)+
  scale_fill_manual(values =  lighter_arm_colors)+
  labs(x = NULL,
       y = "Days to reach threshold of detection",
       color = 'Intervention arm',
       fill = 'Intervention arm') +
  theme_bw(base_size = 14)
