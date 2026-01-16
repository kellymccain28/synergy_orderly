# Plot initial merozoites by intervention arm and by detectability 
plot_initial_merozoites <- function(outputsfolder){
  library(tidyverse)
  library(ggplot2)
  
  path <- 'R:/Kelly/synergy_orderly/src/sim_cohort_generic/outputs/'
  # outputsfolder <- 'outputs_2025-11-26'
  
  # sim_results <- readRDS(paste0(path, outputsfolder, "/sim_results.rds"))
  parasitemia <- readRDS(paste0(path, outputsfolder, "/parasitemia.rds"))#purrr::map_df(sim_results, 'parasitemia')
  # # for now while waiting for output from generic cohort
  # parasitemia <- purrr::map_df(`test_fitted_params_smc_2025-11-24`, 'parasitemia')
  
  initial_values <- parasitemia %>% 
    # Get initial merozoite values 
    filter(time_withinhost == 1) 
  
  ## maybe need to use mero_init???
  VB = 1e6
  detectability_colors <- c('#F4A259', '#5B8E7D')
  lighter <- colorspace::lighten(detectability_colors, amount = 0.3)
  # Adaptation of plot from helper_functions.R make_plots() function  
  ggplot(initial_values %>% filter(mero_init_out!=0)) + 
    geom_boxplot(aes(x = as.factor(det), y = mero_init_out, color = as.factor(det), fill = as.factor(det)), 
                 linewidth = 0.8, alpha = 0.8) + #
    geom_jitter(aes(x = as.factor(det), y = mero_init_out, color = as.factor(det)), alpha = 0.25) + #
    geom_hline(yintercept = 1e-5 * 1e6, color = 'darkred', linetype = 2) +
    facet_wrap(~factor(arm, levels = c('none','rtss','smc','both'))) +
    labs(x = NULL,#'Infection status',
         y = 'Merozoites initating infection') + 
    scale_color_manual(values = detectability_colors) +
    scale_fill_manual(values = lighter) +
    scale_x_discrete(labels = c('0' = 'Cleared', '1' = 'Detectable Case')) +
    theme_bw(base_size = 12) + 
    scale_y_log10(labels = scales::label_log(),
                  breaks = c(1e-9, 1e-7, 1e-5, 1e-3, 0.1, 10, 1000, 100000, 1e7),
                  guide = "axis_logticks"
                  ) +
    theme(legend.position = 'none')
  
  ggsave(paste0(path, outputsfolder,'/initial_merozoites.pdf'), plot = last_plot(),
         height = 5, width = 8)
  
  
  # over all interventions
  ggplot(initial_values) +
    geom_boxplot(aes(x = 1, y = parasites*VB+0001), alpha = 0.7) +
    geom_jitter(aes(x = 1, y = parasites*VB+0.001), alpha = 0.2) + #
    scale_y_log10() +
    facet_wrap(~arm)
  
  initial_values %>% 
    group_by(arm, det) %>%
    summarize(mean_meroinit = mean(mero_init_out),
              mean_pb0 = mean(parasites))
  
  # get percent at 0 in each arm and group 
  initial_values %>%
    group_by(arm) %>%
    mutate(mero_init_out0 = ifelse(mero_init_out == 0, 'start with 0', 'start >0')) %>%
    tabyl(arm, mero_init_out0) %>%
    mutate(perc0 = `start with 0` / (`start with 0`+`start >0`) * 100)
}
