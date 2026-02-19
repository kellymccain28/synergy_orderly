# Plot curve of 1-IRR to understand if there is synergy 
plot_irr <- function(outputsfolder, cohort_folder = 'sim_cohort_generic'){
  # Load packages
  library(zoo)
  library(survival)
  library(survminer)
  library(broom)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  
  source("R:/Kelly/synergy_orderly/shared/format_model_output.R")
  source("R:/Kelly/synergy_orderly/shared/get_incidence.R")
  source("R:/Kelly/synergy_orderly/shared/analyse_model_output.R")
  source("R:/Kelly/synergy_orderly/shared/get_cox_efficacy.R")
  source("R:/Kelly/synergy_orderly/shared/likelihood.R")
  
  path <- paste0('R:/Kelly/synergy_orderly/src/', cohort_folder, '/outputs/')
  # outputsfolder <- 'outputs_2025-12-01_2'
  
  # Using the outputs from monthly_incidence_plot.R
  formatted <- readRDS(paste0(path, outputsfolder, '/formatted_infrecords.rds'))
  inci <- readRDS(paste0(path, outputsfolder, '/incidence.rds'))
  # Using outputs from summarize_IRRs.R
  inci_summary <- readRDS(paste0(path, outputsfolder, '/inci_summary_wide_yearmonth.rds')) %>%
    rename(yearmonth = time_value)
  inci_long <- readRDS(paste0(path, outputsfolder, '/inci_summary_all_yearmonth.rds')) %>%
    filter(metric == 'incidence') %>%
    rename(yearmonth = time_value)

  # Filter for specific date range
  iii_summary <- inci_summary %>% 
    filter(as.Date(yearmonth) > as.Date('2017-06-01') &
             as.Date(yearmonth) < as.Date('2018-01-01'))
  
  metadata_df <- readRDS(paste0(path, outputsfolder, "/metadata_df.rds"))
  base_inputs <- readRDS(paste0(path, outputsfolder, "/base_inputs.rds"))
  params <- readRDS(paste0(path, outputsfolder, "/parameter_grid.rds"))
  if(cohort_folder == 'sim_cohort_generic'){
    if(!is.null(unlist(formatted$smc_dose_days[1]))){
      smc_dates <- as.Date(unlist(formatted$smc_dose_days[1]), origin = '2017-04-01')
    } else if(!is.null(unlist(formatted$smc_dose_days[2]))){
      smc_dates <- as.Date(unlist(formatted$smc_dose_days[2]), origin = '2017-04-01')
    } else if(!is.null(unlist(formatted$smc_dose_days[3]))){
      smc_dates <- as.Date(unlist(formatted$smc_dose_days[3]), origin = '2017-04-01')
    } else if(!is.null(unlist(formatted$smc_dose_days[4]))){
      smc_dates <- as.Date(unlist(formatted$smc_dose_days[4]), origin = '2017-04-01')
    }
  } else if(cohort_folder == 'sim_trial_cohort'){
  smc_dates <- readRDS('R:/Kelly/synergy_orderly/shared/median_smc_dates.rds') %>%
    filter(country == base_inputs$country) %>%
    pull(date)
  }
  smc_lines <- data.frame(
    xintercept = rep(smc_dates,2),
    arm = rep(c('smc', 'both'), each = length(smc_dates)),
    color = '#709176'
  )
  # metadata_df$vaccination_day[1] = 90
  rtss_lines <- data.frame(
    xintercept = as.Date(rep(c(mean(metadata_df$vaccination_day)-60, mean(metadata_df$vaccination_day)-30, mean(metadata_df$vaccination_day), 
                               mean(metadata_df$vaccination_day)[1]+364, mean(metadata_df$vaccination_day)+730),2), origin = '2017-04-01'),
    arm = rep(c('rtss','both'), length(6)),
    color = '#59114D'
  )
  
  
  # First make incidence plot to combine 
  inciall <- ggplot(inci_long) + 
    geom_vline(data = smc_lines, aes(xintercept = xintercept, color = 'SMC delivery'), linetype = 3, linewidth = 0.8)+
    geom_vline(data = rtss_lines, aes(xintercept = xintercept, color = 'RTS,S delivery'), linetype = 2, linewidth = 0.8) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = lower_ci, ymax = upper_ci, fill = arm), 
                alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = median, color = arm), 
              linewidth = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_date(breaks = '3 months',
                 labels = scales::label_date_short()) +
    scale_color_manual(values = c('SMC delivery' = '#709176',
                                  'RTS,S delivery' = '#470024',
                                  'both' = '#E15554', 
                                  'none' = '#E1BC29',
                                  'rtss' = '#3BB273',
                                  'smc' = '#7768AE'),
                       breaks = c('both', 'none', 'rtss', 'smc', 'SMC delivery', 'RTS,S delivery')) +
    scale_fill_manual(values = c('both' = '#E15554', 
                                 'none' = '#E1BC29',
                                 'rtss' = '#3BB273',
                                 'smc' = '#7768AE'),breaks = c('both', 'none', 'rtss', 'smc'),
                      guide = guide_legend(override.aes = list(linetype = 0))) +
    labs(x = 'Date',
         y = 'Incidence per 1000 people',
         color = NULL, fill = NULL) + 
    guides(fill = "none") +
    theme_bw(base_size = 14)
  
  
  # Plot 1: Adding RTSS (filtered)
  ggplot(iii_summary) + 
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = both_smc_q025, ymax = both_smc_q975), 
                fill = '#FFCB77', alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = both_smc_median, color = 'RTS,S added to SMC'), linewidth = 1) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = rtss_none_q025, ymax = rtss_none_q975), 
                fill = '#FE6D73', alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = rtss_none_median, color = 'RTS,S added to none'), linewidth = 1) +
    
    # geom_line(aes(x = as.Date(yearmonth), y = expected_efficacy_median, color = 'Expected efficacy'), linewidth = 1) +
    # geom_line(aes(x = as.Date(yearmonth), y = both_none_median, color = 'Predicted efficacy'), linewidth = 1) +
    ylim(c(0, 1)) + #xlim(c(min(iii_summary$yearmonth), max(iii_summary$yearmonth))) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_date(breaks = '1 month',
                 labels = scales::label_date_short()) +
    scale_color_manual(values = c('RTS,S added to SMC' = '#FFCB77', 
                                  'RTS,S added to none' = '#FE6D73',
                                  'Expected efficacy' = '#6457A6')) +
    labs(x = 'Date',
         y = 'Relative efficacy (1-IRR)',
         color = 'Comparison') + 
    theme_bw(base_size = 14)
  ggsave(paste0(path, outputsfolder,'/irr_addingrtss_filtered.pdf'), plot = last_plot())
  
  # Plot 2: Adding RTSS (non-filtered)
  addrtss <- ggplot(inci_summary) + 
    geom_vline(data = smc_lines, aes(xintercept = xintercept, color = 'SMC delivery'), linetype = 3, linewidth = 0.8)+
    geom_vline(data = rtss_lines, aes(xintercept = xintercept, color = 'RTS,S delivery'), linetype = 2, linewidth = 0.8) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = both_smc_q025, ymax = both_smc_q975), 
                fill = '#FFCB77', alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = both_smc_median, color = 'RTS,S added to SMC'), linewidth = 1) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = rtss_none_q025, ymax = rtss_none_q975), 
                fill = '#FE6D73', alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = rtss_none_median, color = 'RTS,S added to none'), linewidth = 1) +
    
    # geom_line(aes(x = as.Date(yearmonth), y = expected_efficacy_median, color = 'Expected efficacy'), linewidth = 1) +
    
    scale_y_continuous(breaks = c(-0.4,0, 0.25, 0.5, 0.75, 1),
                       limits = c(-0.4, 1)) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_date(breaks = '3 months',
                 labels = scales::label_date_short()) +
    scale_color_manual(values = c('RTS,S added to SMC' = '#FFCB77', 
                                  'RTS,S added to none' = '#FE6D73',
                                  'SMC delivery' = '#709176',
                                  'RTS,S delivery' = '#470024',
                                  'Expected efficacy' = '#6457A6')) +
    labs(x = 'Date',
         y = 'Relative efficacy (1-IRR)',
         color = 'Comparison') + 
    theme_bw(base_size = 14)
  ggsave(paste0(path, outputsfolder,'/irr_addingrtss_nonfiltered.pdf'), plot = addrtss, width = 10, height = 4)
  
  combinedrtss <- plot_grid(addrtss, inciall, nrow = 2,
                            align = 'v')
  ggsave(paste0(path, outputsfolder,'/irr_addingrtss_nonfiltered_andincidence.pdf'), plot = combinedrtss, width = 15, height = 10)
  
  # Filter for SMC comparison
  iii_summary_smc <- inci_summary %>% 
    filter(as.Date(yearmonth) > as.Date('2017-05-01') & as.Date(yearmonth) < as.Date('2018-01-01'))
  
  # Plot 3: Adding SMC (filtered)
  ggplot(iii_summary_smc) + 
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = both_rtss_q025, ymax = both_rtss_q975), 
                fill = '#2ACBCB', alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = both_rtss_median, color = 'SMC added to RTS,S'), linewidth = 1) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = smc_none_q025, ymax = smc_none_q975), 
                fill = '#046865', alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = smc_none_median, color = 'SMC added to none'), linewidth = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    
    # geom_line(aes(x = as.Date(yearmonth), y = expected_efficacy_median, color = 'Expected efficacy'), linewidth = 1) +
    scale_x_date(breaks = '1 month',
                 labels = scales::label_date_short()) +
    scale_color_manual(values = c('SMC added to RTS,S' = '#2ACBCB', 
                                  'SMC added to none' = '#046865',
                                  'Expected efficacy' = '#6457A6')) +
    # ylim(c(-0.5, 1)) +
    labs(x = 'Date',
         y = 'Relative efficacy (1-IRR)',
         color = 'Comparison') + 
    theme_bw(base_size = 14)
  ggsave(paste0(path, outputsfolder,'/irr_addingsmc_filtered.pdf'), plot = last_plot())
  
  # Plot 4: Adding SMC (non-filtered)
  addsmc <- ggplot(inci_summary) + 
    geom_vline(data = smc_lines, aes(xintercept = xintercept, color = 'SMC delivery'), linetype = 2, linewidth = 0.8)+
    geom_vline(data = rtss_lines, aes(xintercept = xintercept, color = 'RTS,S delivery'), linetype = 3, linewidth = 0.8) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = both_rtss_q025, ymax = both_rtss_q975), 
                fill = '#2ACBCB', alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = both_rtss_median, color = 'SMC added to RTS,S'), linewidth = 1) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = smc_none_q025, ymax = smc_none_q975), 
                fill = '#046865', alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = smc_none_median, color = 'SMC added to none'), linewidth = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    
    # geom_line(aes(x = as.Date(yearmonth), y = expected_efficacy_median, color = 'Expected efficacy'), linewidth = 1) +
    scale_x_date(breaks = '3 months',
                 labels = scales::label_date_short()) +
    scale_color_manual(values = c('SMC added to RTS,S' = '#2ACBCB', 
                                  'SMC added to none' = '#046865',
                                  'SMC delivery' = '#709176',
                                  'RTS,S delivery' = '#470024',
                                  'Expected efficacy' = '#6457A6')) +
    # ylim(c(-0.3, 1)) +
    scale_y_continuous(breaks = c(-1, -1.25, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
                       limits = c(-1, 1)) +
    labs(x = 'Date',
         y = 'Relative efficacy (1-IRR)',
         color = 'Comparison') + 
    theme_bw(base_size = 14)
  ggsave(paste0(path, outputsfolder,'/irr_addingsmc_nonfiltered.pdf'), plot = addsmc, width = 10, height = 4)
  
  
  combinedsmc <- plot_grid(addsmc, inciall, nrow = 2,
                          align = 'v')
  ggsave(paste0(path, outputsfolder,'/irr_addingsmc_nonfiltered_andincidence.pdf'), plot = combinedsmc, width = 15, height = 10)
  
  # Print summary statistics
  # cat("Mean of both_smc_irr (filtered):", mean(iii_summary$both_smc_median), "\n")
  # cat("Mean of rtss_none_irr (filtered):", mean(iii_summary$rtss_none_median), "\n")
  # cat("Mean of both_smc_irr (non-filtered):", mean(inci_summary$both_smc_median), "\n")
  # cat("Mean of rtss_none_irr (non-filtered):", mean(inci_summary$rtss_none_median), "\n")
  # 
  # cat("Mean of both_rtss_irr (filtered, Jun 2017-Jan 2018):", mean(iii_summary_smc$both_rtss_median), "\n")
  # cat("Mean of smc_none_irr (filtered, Jun 2017-Jan 2018:", mean(iii_summary_smc$smc_none_median), "\n")
  # cat("Mean of both_rtss_irr (non-filtered):", mean(inci_summary$both_rtss_median), "\n")
  # cat("Mean of smc_none_irr (non-filtered):", mean(inci_summary$smc_none_median), "\n")
  
  
  # Plot of expected vs predicted efficacy of both vs none 
  exppred <- ggplot(inci_summary %>% filter(yearmonth > '2017-05-01')) + 
    geom_vline(data = smc_lines, aes(xintercept = xintercept, color = 'SMC delivery'), linetype = 2, linewidth = 0.8)+
    geom_vline(data = rtss_lines, aes(xintercept = xintercept, color = 'RTS,S delivery'), linetype = 3, linewidth = 0.8) +
    geom_line(aes(x = as.Date(yearmonth), y = expected_efficacy_median, color = 'Expected efficacy'), linewidth = 1) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = expected_efficacy_q025, ymax = expected_efficacy_q975,
                    fill = 'Expected efficacy'), alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = both_none_median, color = 'Model-predicted efficacy'), linewidth = 1) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = both_none_q025, ymax = both_none_q975,
                    fill = 'Model-predicted efficacy'), alpha = 0.3) +
    # ylim(c(-, 1)) + #xlim(c(min(iii_summary$yearmonth), max(iii_summary$yearmonth))) +
    coord_cartesian(ylim = c(-0.15, 1)) +
    scale_y_continuous(breaks=seq(-0.1,1, 0.1),
                       labels=seq(-0.1,1, 0.1)) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_date(breaks = '3 months',
                 labels = scales::label_date_short()) +
    scale_color_manual(values = c('Model-predicted efficacy' = '#59C9A5',
                                  'SMC delivery' = '#709176',
                                  'RTS,S delivery' = '#470024',
                                  'Expected efficacy' = '#6457A6')) +#'#449DD1'
    scale_fill_manual(values = c('Expected efficacy' = '#6457A6',
                                  'Model-predicted efficacy' = '#59C9A5')) +#'#449DD1'
    labs(x = 'Date',
         y = 'Efficacy (1-IRR)',
         color = NULL,
         fill = NULL) + 
    theme_bw(base_size = 14)
  ggsave(paste0(path, outputsfolder,'/predicted_vs_expected_combined_efficacy.pdf'), plot = exppred, width = 10, height = 4)
  
  # Plot of ratio of expected vs predicted efficacy of both vs none 
  ratioplot <-  ggplot(inci_summary %>% filter(yearmonth > '2017-05-01')) + 
    geom_vline(data = smc_lines, aes(xintercept = xintercept, color = 'SMC delivery'), 
               linetype = 2, linewidth = 0.8, alpha = 0.7)+
    geom_vline(data = rtss_lines, aes(xintercept = xintercept, color = 'RTS,S delivery'), 
               linetype = 3, linewidth = 0.8, alpha = 0.7) +
    geom_line(aes(x = as.Date(yearmonth), y = ratio_pred_exp_median), 
              color = '#3E6990', linewidth = 1) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = ratio_pred_exp_q025, ymax = ratio_pred_exp_q975),
                    fill = '#3E6990', alpha = 0.5) +
    # ylim(c(0.2, 1.2)) + #xlim(c(min(iii_summary$yearmonth), max(iii_summary$yearmonth))) +
     # coord_cartesian(ylim = c(0.8,1.2)) + 
    # scale_y_continuous(breaks = seq(-0.4,1.2,0.1),
    #                    labels = seq(-0.4,1.2,0.1))+
    geom_hline(yintercept = 1, linetype = 2) +
    scale_x_date(breaks = '3 months',
                 labels = scales::label_date_short()) +
    scale_color_manual(values = c('SMC delivery' = '#709176',
                                  'RTS,S delivery' = '#470024')) +#'#449DD1'
    # scale_fill_manual(values = c('Expected efficacy' = '#6457A6',
    #                              'Model-predicted efficacy' = '#59C9A5')) +#'#449DD1'
    labs(x = 'Date',
         y = 'Ratio of model-predicted to expected efficacy',
         color = NULL,
         fill = NULL) + 
    theme_bw(base_size = 14)
   
   combinedratio <- plot_grid(exppred, ratioplot + theme(legend.position = 'none'), nrow = 2,
                             align = 'v')
   
   # mean(inci_summary$ratio_pred_exp_median)
   # mean(inci_summary$ratio_pred_exp_q025)
   # mean(inci_summary$ratio_pred_exp_q975)
  ggsave(paste0(path, outputsfolder,'/predicted_vs_expected_combined_ratio.pdf'), plot = ratioplot, width = 10, height = 4)
  ggsave(paste0(path, outputsfolder,'/predicted_vs_expected_combined_and_ratio.pdf'), plot = combinedratio, width = 12, height = 9)
  
  
  # Plot of expected vs predicted cases averted per 1000 of both vs none 
  exppred_inci <- ggplot(inci_summary %>% filter(yearmonth > '2017-03-01')) + 
    geom_vline(data = smc_lines, aes(xintercept = xintercept, color = 'SMC delivery'), linetype = 2, linewidth = 0.8)+
    geom_vline(data = rtss_lines, aes(xintercept = xintercept, color = 'RTS,S delivery'), linetype = 3, linewidth = 0.8) +
    geom_line(aes(x = as.Date(yearmonth), y = inci_averted_expected_median, color = 'Expected cases averted\nper 1000 person months'), linewidth = 1) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = inci_averted_expected_q025, ymax = inci_averted_expected_q975,
                    fill = 'Expected cases averted\nper 1000 person months'), alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = inci_averted_model_median, color = 'Model-predicted cases averted\nper 1000 person months'), 
              linewidth = 1) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = inci_averted_model_q025, ymax = inci_averted_model_q975,
                    fill = 'Model-predicted cases averted\nper 1000 person months'), alpha = 0.3) +
    # ylim(c(-, 1)) + #xlim(c(min(iii_summary$yearmonth), max(iii_summary$yearmonth))) +
    # coord_cartesian(ylim = c(-0.15, 1)) +
    # scale_y_continuous(breaks=seq(-0.1,1, 0.1),
    #                    labels=seq(-0.1,1, 0.1)) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_date(breaks = '3 months',
                 labels = scales::label_date_short()) +
    scale_color_manual(values = c('Model-predicted cases averted\nper 1000 person months' = '#59C9A5',
                                  'SMC delivery' = '#709176',
                                  'RTS,S delivery' = '#470024',
                                  'Expected cases averted\nper 1000 person months' = '#6457A6')) +#'#449DD1'
    scale_fill_manual(values = c('Expected cases averted\nper 1000 person months' = '#6457A6',
                                 'Model-predicted cases averted\nper 1000 person months' = '#59C9A5')) +#'#449DD1'
    labs(x = 'Date',
         y = 'Cases averted per 1000 person months',
         color = NULL,
         fill = NULL) + 
    theme_bw(base_size = 14)
  ggsave(paste0(path, outputsfolder,'/predicted_vs_expected_combined_incidenceaverted.pdf'), plot = exppred_inci, width = 10, height = 4)
  
  # Plot of difference of expected vs predicted cases averted per 1000 of both vs none 
  diffplot <-  ggplot(inci_summary %>% filter(yearmonth > '2017-05-01')) + 
    geom_vline(data = smc_lines, aes(xintercept = xintercept, color = 'SMC delivery'), 
               linetype = 2, linewidth = 0.8, alpha = 0.7)+
    geom_vline(data = rtss_lines, aes(xintercept = xintercept, color = 'RTS,S delivery'), 
               linetype = 3, linewidth = 0.8, alpha = 0.7) +
    geom_line(aes(x = as.Date(yearmonth), y = difference_inci_averted_pred_exp_median), 
              color = '#3E6990', linewidth = 1) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = difference_inci_averted_pred_exp_q025, 
                    ymax = difference_inci_averted_pred_exp_q975),
                fill = '#3E6990', alpha = 0.5) +
    geom_hline(yintercept = 1, linetype = 2) +
    scale_x_date(breaks = '3 months',
                 labels = scales::label_date_short()) +
    scale_color_manual(values = c('SMC delivery' = '#709176',
                                  'RTS,S delivery' = '#470024')) +#'#449DD1'
    # scale_fill_manual(values = c('Expected efficacy' = '#6457A6',
    #                              'Model-predicted efficacy' = '#59C9A5')) +#'#449DD1'
    labs(x = 'Date',
         y = 'Difference in model-predicted cases averted per 1000\nversus expected cases averted per 1000',
         color = NULL,
         fill = NULL) + 
    theme_bw(base_size = 14)
  
  combineddifference <- plot_grid(exppred_inci, diffplot + theme(legend.position = 'none'), nrow = 2,
                             align = 'v')
  ggsave(paste0(path, outputsfolder,'/predicted_vs_expected_combined_incidence.pdf'), plot = diffplot, width = 10, height = 4)
  ggsave(paste0(path, outputsfolder,'/predicted_vs_expected_combined_incidence_and_difference.pdf'), plot = combineddifference, width = 12, height = 9)
  
  
  # Find percentage of cases that fall in July to November each year 
  cases_in_season <- inci %>% 
    mutate(in_season = case_when(
      (yearmonth >= '2017-06-01' & yearmonth <= '2017-11-30') ~ 1, 
      (yearmonth >= '2018-06-01' & yearmonth <= '2018-11-30') ~ 1,
      (yearmonth >= '2019-06-01' & yearmonth <= '2019-11-30') ~ 1,
      TRUE ~ 0)) %>%
    group_by(in_season) %>%
    summarize(n_cases = sum(n_cases)) %>%
    mutate(total = sum(n_cases),
           p_in_season = n_cases / total * 100)
  cases_in_season
}