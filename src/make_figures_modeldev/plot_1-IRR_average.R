# plot overall 1-IRR for each year 

plot_irr_average <- function(outputsfolder, 
                             agg_unit = 'year',
                             cohort_folder = 'sim_cohort_generic'){
  # Load packages
  library(zoo)
  library(survival)
  library(survminer)
  library(broom)
  library(ggplot2)
  library(tidyverse)
  
  source("R:/Kelly/synergy_orderly/shared/format_model_output.R")
  source("R:/Kelly/synergy_orderly/shared/get_incidence.R")
  source("R:/Kelly/synergy_orderly/shared/analyse_model_output.R")
  source("R:/Kelly/synergy_orderly/shared/get_cox_efficacy.R")
  source("R:/Kelly/synergy_orderly/shared/likelihood.R")
  
  path <- paste0('R:/Kelly/synergy_orderly/src/', cohort_folder, '/outputs/')
  
  # Using the outputs from monthly_incidence_plot.R
  # inci <- readRDS(paste0(path, outputsfolder, '/incidence.rds'))
  
  # Using outputs from summarize_IRRs.R
  inci_summary <- readRDS(paste0(path, outputsfolder, '/inci_summary_all_', agg_unit, '.rds'))
  
  # get the values for the text
  inci_summary %>%
    filter(time_value == 'overall' & metric == 'efficacy')  %>%
    saveRDS(paste0(path, outputsfolder, '/summary_efficacy.rds'))
  
  # Plotting
  colors <- RColorBrewer::brewer.pal(9, 'BuPu')
  colors <- c(colors[4:length(colors)], 'black')
  
  inci_summary <- inci_summary %>%
    mutate(shape_var = case_when(
      metric == 'efficacy' & grepl('Expected', comparison) ~ 'Expected',
      metric == 'efficacy' ~ 'Model-predicted',
      TRUE ~ NA)) %>%
    # rbind(expected %>% 
    #         mutate(comparison = 'Expected: both vs none',
    #                shape_var = 'Expected')) %>%
    mutate(comparison = factor(comparison, levels = c("Expected both vs none", 'both vs none',
                                                      'both vs rtss','both vs smc',
                                                      'rtss vs none','smc vs none',
                                                      'rtss vs smc','smc vs rtss')))
  # plot
  ggplot(inci_summary %>% filter(metric == 'efficacy'), aes(x = comparison, y = median, 
                                                            color = time_value, group = time_value, shape = shape_var)) +
    # model estimated
    geom_point(position = position_dodge(width = 0.5), size = 1) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                  width = 0.2, position = position_dodge(width = 0.5)) +
    scale_shape_manual(values = c('Expected' = 7, 'Model-predicted' = 16)) + 
    # scale_linetype_manual(values = c('Expected: both vs none' = 2)) +
    # scale_color_brewer(palette = 'BuPu') +
    scale_color_manual(values = colors) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
    labs(
      x = "Arm Comparison",
      y = "Relative efficacy (1-IRR)",
      # title = "Median IRR with 95% CI by Intervention Comparison",
      shape = NULL, linetype = NULL,
      color = if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Study half-year'
    ) +
    scale_y_continuous(breaks = seq(-1,1,0.2)) + 
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  ggsave(paste0(path, outputsfolder,'/irr_average_by', agg_unit, '.pdf'), plot = last_plot(), width = 12, height = 6)
    
    ggplot(inci_summary %>% filter(!is.na(shape_var) & grepl('oth vs none', comparison)), 
           aes(x = time_value, y = median, color = shape_var, group = shape_var)) +
      # model estimated
      geom_point(aes(shape = shape_var), position = position_dodge(width = 0.5)) +
      geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                    width = 0.2, position = position_dodge(width = 0.5)) +
      scale_shape_manual(values = c('Expected' = 7, 'Model-predicted' = 16))+#, guide = NULL) +
      scale_linetype_manual(values = c('Expected' = 2)) +
      scale_color_brewer(palette = 'Dark2') +
      labs(
        x = if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Half-year',
        y = "Relative efficacy (1-IRR)",
        shape = NULL, linetype = NULL,
        color = NULL
      ) +
      # scale_y_continuous(breaks = seq(0,1,0.2)) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0(path, outputsfolder,'/bothvsnone_irr_', agg_unit, '.pdf'), plot = last_plot(), width = 6, height = 4)
    
    
    # Plot of the ratio of model-predicted to expected by aggregation unit 
    ggplot(inci_summary %>% filter(metric == 'ratio pred to exp'), 
           aes(x = time_value, y = median, 
               color = time_value, group = time_value)) +
      # model estimated
      geom_point(size = 1) +
      geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                    width = 0.2, position = position_dodge(width = 0.5)) +
      scale_color_manual(values = colors) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
      labs(
        x = if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Half-year',
        y = "Ratio of model-predicted to expected efficacy\n of combination vs no intervention",
        # title = "Median IRR with 95% CI by Intervention Comparison",
        shape = NULL, linetype = NULL,
        color = NULL#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Study half-year'
      ) +
      # scale_y_continuous(breaks = seq(-1,1,0.2)) + 
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = 'none')
    
    ggsave(paste0(path, outputsfolder,'/ratio_average_by', agg_unit, '.pdf'), plot = last_plot(), width = 8, height = 6)
    inci_summary %>% filter(metric == 'ratio pred to exp')

   
    
}
