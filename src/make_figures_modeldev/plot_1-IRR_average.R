# plot overall 1-IRR for each year 

plot_irr_average <- function(outputsfolder, # in thesis, have used may 1 and aug 26, scen 5, outputs_2026-02-18_3
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
    filter(time_value == 'Overall' & (metric == 'efficacy' | 
                                        metric == 'cases_averted_model'))  %>%
    mutate(median = median * 100,
           lower_ci = lower_ci * 100,
           upper_ci = upper_ci * 100) %>%
    dplyr::select(comparison, median, lower_ci, upper_ci) %>%
    saveRDS(paste0(path, outputsfolder, '/summary_efficacy.rds'))
  
  # Plotting
  colors <- RColorBrewer::brewer.pal(9, 'BuPu')
  colors <- c(colors[4:length(colors)], 'black')
  
  inci_summary <- inci_summary %>%
    mutate(shape_var = case_when(
      metric == 'efficacy' & grepl('Expected', comparison) ~ 'Expected',
      metric == 'efficacy' ~ 'Model-predicted',
      TRUE ~ NA)) %>%
    mutate(comparison = factor(comparison, levels = c("Expected both vs none", 'both vs none',
                                                      'both vs rtss','both vs smc',
                                                      'rtss vs none','smc vs none',
                                                      'rtss vs smc','smc vs rtss')))
  # plot
  # ggplot(inci_summary %>% filter(metric == 'efficacy'), aes(x = comparison, y = median, 
  #                                                           color = time_value, group = time_value, shape = shape_var)) +
  #   # model estimated
  #   geom_point(position = position_dodge(width = 0.5), size = 1) +
  #   geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
  #                 width = 0.2, position = position_dodge(width = 0.5)) +
  #   scale_shape_manual(values = c('Expected' = 7, 'Model-predicted' = 16)) + 
  #   # scale_linetype_manual(values = c('Expected: both vs none' = 2)) +
  #   # scale_color_brewer(palette = 'BuPu') +
  #   scale_color_manual(values = colors) +
  #   geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
  #   labs(
  #     x = "Arm Comparison",
  #     y = "Relative efficacy (1-IRR)",
  #     # title = "Median IRR with 95% CI by Intervention Comparison",
  #     shape = NULL, linetype = NULL,
  #     color = if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Study half-year'
  #   ) +
  #   scale_y_continuous(breaks = seq(-1,1,0.2)) + 
  #   theme_minimal(base_size = 14) +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  #   
  # ggsave(paste0(path, outputsfolder,'/irr_average_by', agg_unit, '.pdf'), plot = last_plot(), width = 12, height = 6)
    
  
  # irr overall 
  # plot
  ggplot(inci_summary %>% filter(metric == 'efficacy' & time_value == 'Overall'), 
         aes(x = comparison, y = median, 
             shape = shape_var, color = shape_var)) +
    # model estimated
    geom_point(position = position_dodge(width = 0.5), size = 1) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                  width = 0.2, position = position_dodge(width = 0.5)) +
    scale_shape_manual(values = c('Expected' = 7, 'Model-predicted' = 16)) + 
    scale_color_manual(values = c('Expected' = "#8C6BB1", 'Model-predicted' ="#810F7C")) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
    labs(
      x = "Arm Comparison",
      y = "Relative efficacy (1-IRR)",
      shape = NULL, linetype = NULL,
      color = NULL#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Study half-year'
    ) +
    scale_y_continuous(breaks = seq(-1,1,0.2)) + 
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(paste0(path, outputsfolder,'/irr_average_overall.pdf'), plot = last_plot(), width = 8, height = 6)
  
  
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
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0(path, outputsfolder,'/bothvsnone_irr_', agg_unit, '.pdf'), plot = last_plot(), width = 8, height = 6)
    
    
    # Plot of the difference of model-predicted to expected by aggregation unit 
    ggplot(inci_summary %>% filter(metric == 'difference_inci_averted_pred_exp'), 
           aes(x = time_value, y = median, 
               color = time_value, group = time_value)) +
      # model estimated
      geom_point(size = 1) +
      geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                    width = 0.2, position = position_dodge(width = 0.5)) +
      scale_color_manual(values = colors) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
      # Add annotations for synergistic/antagonistic
      annotate("text", x = -Inf, y = max(inci_summary$upper_ci[inci_summary$metric == 'difference_inci_averted_pred_exp'], na.rm = TRUE) * 0.05, 
               label = "Synergistic", color = "green4", hjust = -0.1, vjust = 1, size = 4, fontface = "italic") +
      annotate("text", x = -Inf, y = min(inci_summary$lower_ci[inci_summary$metric == 'difference_inci_averted_pred_exp'], na.rm = TRUE) * 0.05, 
               label = "Antagonistic", color = "red3", hjust = -0.1, vjust = 1.5, size = 4, fontface = "italic") +
      labs(
        x = NULL,#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Half-year',
        y = "Difference in model-predicted versus expected\ncases averted per 1000 person-months of\ncombination vs no intervention",
        shape = NULL, linetype = NULL,
        color = NULL#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Study half-year'
      ) +
      # scale_y_continuous(breaks = seq(-1,1,0.2)) + 
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = 'none')
    
    ggsave(paste0(path, outputsfolder,'/difference_inci_average_by', agg_unit, '.pdf'), plot = last_plot(), width = 8, height = 6)
    inci_summary %>% filter(metric == 'difference_inci_averted_pred_exp')
    
    # print percentage difference
    inci_comparison <- inci_summary %>%
      filter(metric == 'inci_averted_model' | metric == 'inci_averted_expected' | 
               metric == 'difference_inci_averted_pred_exp' | metric == 'pct_diff_inci_averted') %>%
      dplyr::select(time_value, time_value_num, metric, median) %>%
      pivot_wider(id_cols = c(time_value, time_value_num), 
                  names_from = metric, 
                  values_from = median) %>%
      mutate(
        abs_diff = inci_averted_model - inci_averted_expected
      )
    
    inci_summary %>%
      filter(metric == 'difference_inci_averted_pred_exp')
    
    saveRDS(inci_comparison, paste0(path, outputsfolder,'/inci_percent_averted.rds'))
    
}
