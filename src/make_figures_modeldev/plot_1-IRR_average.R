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
  
  # if(agg_unit == 'year'){
  #   # Get annual and overall incidence
  # inci_annual <- inci %>%
  #   filter(!is.na(date)) %>%
  #   mutate(studyyear = case_when(date < '2018-04-01' ~ 1,
  #                                date < '2019-04-01' ~ 2,
  #                                date < '2020-04-01' ~ 3)) %>%
  #   group_by(studyyear, arm, sim_id) %>%
  #   summarize(person_months = sum(person_months),
  #             n_cases = sum(n_cases)) %>%
  #   mutate(incidence_per_1000pm = n_cases / person_months * 1000,
  #          studyyear = as.character(studyyear))
  # } else if (agg_unit == 'halfyear'){
  #   inci_annual <- inci %>%
  #     filter(!is.na(date)) %>%
  #     mutate(studyyear = case_when(date < '2017-10-01' ~ 1,
  #                             date < '2018-04-01' ~ 2,
  #                             date < '2018-10-01' ~ 3,
  #                             date < '2019-04-01' ~ 4,
  #                             date < '2019-10-01' ~ 5,
  #                             date < '2020-04-01' ~ 6)) %>%
  #     group_by(studyyear, arm, sim_id) %>%
  #     summarize(person_months = sum(person_months),
  #               n_cases = sum(n_cases)) %>%
  #     mutate(incidence_per_1000pm = n_cases / person_months * 1000,
  #            studyyear = as.character(studyyear))
  # }
  # 
  # inci_overall <- inci %>%
  #   group_by(arm, sim_id) %>%
  #   summarize(person_months = sum(person_months),
  #             n_cases = sum(n_cases)) %>%
  #   mutate(incidence_per_1000pm = n_cases / person_months * 1000,
  #          studyyear = 'overall')
  # 
  # inci_averages <- rbind(inci_annual, inci_overall)
  # 
  # # Pivot wider
  # inci_wide <- inci_averages %>%
  #   split(.$sim_id) %>%
  #   map_dfr(~ .x %>%
  #             select(arm, studyyear,
  #                    person_months, incidence_per_1000pm) %>%
  #             pivot_wider(
  #               names_from = arm,
  #               values_from = c(person_months, incidence_per_1000pm),
  #               id_cols = c(studyyear)
  #             ),
  #           .id = "sim_id")
  # 
  # # Get IRRs
  # inci_wide <- inci_wide %>%
  #   mutate(rtss_none_irr = (incidence_per_1000pm_rtss / incidence_per_1000pm_none),
  #          smc_none_irr = (incidence_per_1000pm_smc / incidence_per_1000pm_none),
  #          both_none_irr = (incidence_per_1000pm_both / incidence_per_1000pm_none),
  #          rtss_smc_irr = (incidence_per_1000pm_rtss / incidence_per_1000pm_smc),
  #          both_smc_irr = (incidence_per_1000pm_both / incidence_per_1000pm_smc),
  #          both_rtss_irr = (incidence_per_1000pm_both / incidence_per_1000pm_rtss),
  #          smc_rtss_irr = (incidence_per_1000pm_smc / incidence_per_1000pm_rtss)  )%>%
  #   mutate(expected_efficacy = 1 - (rtss_none_irr * smc_none_irr))
  # 
  # # Calculate median and IQR for each aggregation unit
  # inci_summary <- inci_wide %>%
  #   group_by(studyyear) %>%
  #   summarise(
  #     # Both vs SMC
  #     both_smc_median = median(both_smc_irr, na.rm = TRUE),
  #     both_smc_q25 = quantile(both_smc_irr, 0.25, na.rm = TRUE),
  #     both_smc_q75 = quantile(both_smc_irr, 0.75, na.rm = TRUE),
  #     # RTSS vs none
  #     rtss_none_median = median(rtss_none_irr, na.rm = TRUE),
  #     rtss_none_q25 = quantile(rtss_none_irr, 0.25, na.rm = TRUE),
  #     rtss_none_q75 = quantile(rtss_none_irr, 0.75, na.rm = TRUE),
  #     # Both vs RTSS
  #     both_rtss_median = median(both_rtss_irr, na.rm = TRUE),
  #     both_rtss_q25 = quantile(both_rtss_irr, 0.25, na.rm = TRUE),
  #     both_rtss_q75 = quantile(both_rtss_irr, 0.75, na.rm = TRUE),
  #     # SMC vs none
  #     smc_none_median = median(smc_none_irr, na.rm = TRUE),
  #     smc_none_q25 = quantile(smc_none_irr, 0.25, na.rm = TRUE),
  #     smc_none_q75 = quantile(smc_none_irr, 0.75, na.rm = TRUE),
  #     .groups = 'drop'
  #   )
  # 
  # # Reshape the data to long format for plotting
  # inci_long <- inci_wide %>%
  #   pivot_longer(
  #     cols = ends_with("_irr"),
  #     names_to = "comparison",
  #     values_to = "irr"
  #   ) %>%
  #   mutate(
  #     # Clean up the comparison names for better labels
  #     comparison = gsub("_irr$", "", comparison),
  #     comparison = gsub("_", " vs ", comparison)
  #   )
  # 
  # # Calculate median and 95% CI for each comparison
  # inci_summary <- inci_long %>%
  #   group_by(comparison, studyyear) %>%
  #   summarise(
  #     median = median(1 - irr, na.rm = TRUE),
  #     lower_ci = quantile(1 - irr, 0.025, na.rm = TRUE),
  #     upper_ci = quantile(1 - irr, 0.975, na.rm = TRUE)
  #   )
  # 
  # # Calculate the expected value of IRR
  # rtssnone <- inci_summary[inci_summary$comparison == 'rtss vs none', c('median', "lower_ci", "upper_ci")]
  # smcnone <- inci_summary[inci_summary$comparison == 'smc vs none', c('median', "lower_ci", "upper_ci")]#$VE
  # years <- inci_summary[inci_summary$comparison == 'rtss vs none',]$studyyear
  # 
  # expected <- 1 - ((1-rtssnone) * (1-smcnone)) #rtssnone * (1 - smcnone) + smcnone -- second equation is the one from Sherrard-Smith
  # expected$studyyear <- inci_summary[inci_summary$comparison == 'rtss vs none',]$studyyear
  
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
  ggplot(inci_summary %>% filter(metric == 'efficacy'), aes(x = comparison, y = median, color = time_value, group = time_value, shape = shape_var)) +
    # model estimated
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                  width = 0.2, position = position_dodge(width = 0.5)) +
    # expected 
    # geom_point(data = expected, 
    #            aes(x = 'both vs none', y = median, color = studyyear, group = studyyear, shape = 'Expected'),
    #            position = position_dodge(width = 0.8)) +
    # geom_errorbar(data = expected, 
    #               aes(x = 'both vs none', y = median, color = studyyear, group = studyyear, ymin = lower_ci, ymax = upper_ci,
    #                   linetype = 'Expected'),
    #               width = 0.2, position = position_dodge(width = 0.8)) +
    
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
      color = 'Study year'
    ) +
    scale_y_continuous(breaks = seq(-1,1,0.2)) + 
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(paste0(path, outputsfolder,'/irr_average_by', agg_unit, '.pdf'), plot = last_plot(), width = 12, height = 6)
    
    
    # Synergy plot
    expected <- inci_summary %>%
      filter(comparison == 'Expected both vs none') %>%
      mutate(median_exp = median,
             lower_ci_exp = lower_ci, 
             upper_ci_exp = upper_ci)
    
    inci_summary2 <- inci_summary %>%
      separate(comparison, 
               into = c("comparator", "reference"), 
               sep = " vs ",
               remove = FALSE) %>%
      left_join(expected) %>%
      mutate(synergistic = ifelse(median > median_exp, 1, 0))
    
    # ggplot(inci_summary2)+
    #   geom_tile(aes(x = comparator, y = reference, fill = median)) +
    #   scale_fill_distiller(palette = 'RdYlGn', direction = 1) +
    #   facet_wrap(~time_value)
    
    # syn <- inci_summary2 %>%
    #   filter(comparison == 'both vs none') %>%
    #   mutate(expected = 'Model-predicted') %>%
    #   rbind(expected %>% mutate(expected = 'Expected'))
    
    # ggplot(syn %>% filter(!is.na(shape_var)), aes(x = time_value, y = median, color = expected, group = expected)) +
    #   # model estimated
    #   geom_point(aes(shape = expected), position = position_dodge(width = 0.5)) +
    #   geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
    #                 width = 0.2, position = position_dodge(width = 0.5)) +
    #   scale_shape_manual(values = c('Expected' = 7, 'Model-predicted' = 16))+#, guide = NULL) +
    #   scale_linetype_manual(values = c('Expected' = 2)) +
    #   scale_color_brewer(palette = 'Dark2') +
    #   # scale_color_manual(values = colors) +
    #   geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
    #   labs(
    #     x = if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Half-year',
    #     y = "Relative efficacy (1-IRR)",
    #     # title = "Median IRR with 95% CI by Intervention Comparison",
    #     shape = NULL, linetype = NULL,
    #     color = NULL
    #   ) +
    #   scale_y_continuous(breaks = seq(0,1,0.2)) +
    #   theme_minimal(base_size = 12) +
    #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggplot(inci_summary %>% filter(!is.na(shape_var) & grepl('oth vs none', comparison)), aes(x = time_value, y = median, color = shape_var, group = shape_var)) +
      # model estimated
      geom_point(aes(shape = shape_var), position = position_dodge(width = 0.5)) +
      geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                    width = 0.2, position = position_dodge(width = 0.5)) +
      scale_shape_manual(values = c('Expected' = 7, 'Model-predicted' = 16))+#, guide = NULL) +
      scale_linetype_manual(values = c('Expected' = 2)) +
      scale_color_brewer(palette = 'Dark2') +
      # scale_color_manual(values = colors) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
      labs(
        x = if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Half-year',
        y = "Relative efficacy (1-IRR)",
        # title = "Median IRR with 95% CI by Intervention Comparison",
        shape = NULL, linetype = NULL,
        color = NULL
      ) +
      scale_y_continuous(breaks = seq(0,1,0.2)) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0(path, outputsfolder,'/bothvsnone_irr_', agg_unit, '.pdf'), plot = last_plot(), width = 6, height = 4)
    
    
    # Plot of the ratio of model-predicted to expected by aggregation unit 
    

    # # Plot of only values overall 
    # ggplot(inci_summary %>% filter(studyyear == 'overall'), 
    #        aes(x = comparison, y = median)) +
    #   # model estimated
    #   geom_point(aes(color = reference, group = reference), 
    #              position = position_dodge(width = 0.5)) +
    #   geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci, color = reference, group = reference),
    #                 width = 0.2, position = position_dodge(width = 0.5)) +
    #   # expected 
    #   geom_point(data = expected %>% filter(studyyear == 'overall'), 
    #              aes(x = 'both vs none', y = median, shape = 'Expected'),
    #              color = colorspace::lighten('#D96C06', amount = 0.4),
    #              position = position_dodge(width = 2)) +
    #   geom_errorbar(data = expected %>% filter(studyyear == 'overall'), 
    #                 aes(x = 'both vs none', y = median, ymin = lower_ci, ymax = upper_ci,
    #                     linetype = 'Expected'), color = colorspace::lighten('#D96C06', amount = 0.4),
    #                 width = 0.2, position = position_dodge(width = 2)) +
    #   scale_shape_manual(values = c('Expected' = 7)) + 
    #   scale_linetype_manual(values = c('Expected' = 2)) +
    #   # scale_color_brewer(palette = 'BuPu') +
    #   scale_color_brewer(palette = 'Dark2') +
    #   geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
    #   labs(
    #     x = "Arm Comparison",
    #     y = "Relative efficacy (1-IRR)",
    #     # title = "Median IRR with 95% CI by Intervention Comparison",
    #     shape = NULL, linetype = NULL,
    #     color = NULL
    #   ) +
    #   scale_y_continuous(breaks = seq(-1,1,0.2)) + 
    #   theme_minimal(base_size = 12) +
    #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # ggsave(paste0(path, outputsfolder,'/irr_average_overall.pdf'), plot = last_plot(), width = 12, height = 6)
    
    
}
