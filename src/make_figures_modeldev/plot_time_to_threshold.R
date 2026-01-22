# Script to plot the time it takes for each infection to reach the threshold, by arm 
plot_time_to_threshold <- function(outputsfolder){
  # Load packages
  library(ggplot2)
  library(colorspace)
  library(janitor)
  
  t_liverstage = 8
  arm_colors <- c('#E1BC29',
                  '#7768AE',
                  '#3BB273',
                  '#E15554')
  lighter_arm_colors <- colorspace::lighten(arm_colors, amount = 0.3)
  
  path <- 'R:/Kelly/synergy_orderly/src/sim_cohort_generic/outputs/'

  # Using the outputs from monthly_incidence_plot.R
  formatted <- readRDS(paste0(path, outputsfolder, '/formatted_infrecords.rds')) %>%
    filter(detectable==1 ) %>%
    mutate(arm = factor(arm, levels = c('none','smc','rtss','both'))) %>%
    mutate(t_toreach_threshold = t_toreach_threshold + t_liverstage) # to account for liver stage time 
  
  distr <- formatted %>% filter(sim_id == 'parameter_set_1_generic_0.9') %>% tabyl(arm)
  
  ggplot(formatted %>% filter(sim_id == 'parameter_set_1_generic_0.9')) + 
    geom_boxplot(aes(x = arm, y = t_toreach_threshold, color = arm, fill = arm), 
                 linewidth = 0.8, outlier.alpha = 0.1) + 
    geom_text(data = distr, 
              aes(x = arm, y = 12, label = paste0('n = ', round(n, 1))),
                  size = 3.5, fontface = 'bold') +
    # geom_jitter(aes(x = arm, y = t_toreach_threshold, color = arm), alpha = 0.2) + 
    scale_y_log10(breaks = c(10, 30, 50, 100, 200, 300),
                  labels = c(10, 30, 50, 100, 200, 300),
                  limits = c(10,NA),
                  guide = "axis_logticks") + 
    scale_color_manual(values =  arm_colors)+
    scale_fill_manual(values =  lighter_arm_colors)+
    labs(x = NULL,
         y = "Days to reach clinical threshold",
         color = 'Intervention arm',
         fill = 'Intervention arm') +
    theme_bw(base_size = 14) +
    theme(legend.position = 'none') 
  ggsave('R:/Kelly/synergy_orderly/figures/time_to_threshold.pdf', plot = last_plot(), height = 4, width = 6)
  ggsave(paste0(path, outputsfolder,'/time_to_threshold.pdf'), plot = last_plot(), height = 4, width = 6)
  
  formatted %>%
    group_by(arm) %>%
    summarize(min_t_threshold = min(t_toreach_threshold),
              max_t_threshold = max(t_toreach_threshold),
              median_t_threshold = median(t_toreach_threshold),
              mean_t_threshold = mean(t_toreach_threshold),
              lower_t_threshold = quantile(t_toreach_threshold, 0.025),
              upper_t_threshold = quantile(t_toreach_threshold, 0.975))
}
# A tibble: 4 × 5
# arm   median_t_threshold mean_t_threshold lower_t_threshold upper_t_threshold
# <fct>              <dbl>            <dbl>             <dbl>             <dbl>
# 1 none                  26             26.4                22                32
# 2 smc                   36             59.7                24               160
# 3 rtss                  28             28.0                24                34
# 4 both                  38             56.9                24               152