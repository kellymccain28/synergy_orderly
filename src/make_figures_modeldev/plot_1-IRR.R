# Plot curve of 1-IRR to understand if there is synergy 
plot_irr <- function(outputsfolder){
  # Load packages
  library(zoo)
  library(survival)
  library(survminer)
  library(broom)
  library(ggplot2)
  
  source("R:/Kelly/synergy_orderly/shared/format_model_output.R")
  source("R:/Kelly/synergy_orderly/shared/get_incidence.R")
  source("R:/Kelly/synergy_orderly/shared/analyse_model_output.R")
  source("R:/Kelly/synergy_orderly/shared/get_cox_efficacy.R")
  source("R:/Kelly/synergy_orderly/shared/likelihood.R")
  
  path <- 'R:/Kelly/synergy_orderly/src/sim_cohort_generic/outputs/'
  # outputsfolder <- 'outputs_2025-12-01_2'
  
  # Using the outputs from monthly_incidence_plot.R
  formatted <- readRDS(paste0(path, outputsfolder, '/formatted_infrecords.rds'))
  inci <- readRDS(paste0(path, outputsfolder, '/incidence.rds'))
  
  # First need to get IRR values for each combination 
  
  # Pivot wider 
  inci_wide <- inci %>%
    split(.$sim_id) %>%
    map_dfr(~ .x %>%
              select(arm, year, month, yearmonth,
                     person_months, rate, incidence_per_1000pm) %>%
              pivot_wider(
                names_from = arm,
                values_from = c(person_months, rate, incidence_per_1000pm),
                id_cols = c(year, month, yearmonth)
              ),
            .id = "sim_id")
  
  
  # Get IRRs
  inci_wide <- inci_wide %>%
    mutate(rtss_none_irr = 1 - (incidence_per_1000pm_rtss / incidence_per_1000pm_none),
           smc_none_irr = 1 - (incidence_per_1000pm_smc / incidence_per_1000pm_none),
           both_none_irr = 1 - (incidence_per_1000pm_both / incidence_per_1000pm_none),
           rtss_smc_irr = 1 - (incidence_per_1000pm_rtss / incidence_per_1000pm_smc),
           both_smc_irr = 1 - (incidence_per_1000pm_both / incidence_per_1000pm_smc),
           both_rtss_irr = 1 - (incidence_per_1000pm_both / incidence_per_1000pm_rtss),
           smc_rtss_irr = 1 - (incidence_per_1000pm_both / incidence_per_1000pm_rtss)  )
  
  # Calculate median and IQR for each month
  inci_summary <- inci_wide %>%
    group_by(yearmonth) %>%
    summarise(
      # Both vs SMC
      both_smc_median = median(both_smc_irr, na.rm = TRUE),
      both_smc_q25 = quantile(both_smc_irr, 0.25, na.rm = TRUE),
      both_smc_q75 = quantile(both_smc_irr, 0.75, na.rm = TRUE),
      # RTSS vs none
      rtss_none_median = median(rtss_none_irr, na.rm = TRUE),
      rtss_none_q25 = quantile(rtss_none_irr, 0.25, na.rm = TRUE),
      rtss_none_q75 = quantile(rtss_none_irr, 0.75, na.rm = TRUE),
      # Both vs RTSS
      both_rtss_median = median(both_rtss_irr, na.rm = TRUE),
      both_rtss_q25 = quantile(both_rtss_irr, 0.25, na.rm = TRUE),
      both_rtss_q75 = quantile(both_rtss_irr, 0.75, na.rm = TRUE),
      # SMC vs none
      smc_none_median = median(smc_none_irr, na.rm = TRUE),
      smc_none_q25 = quantile(smc_none_irr, 0.25, na.rm = TRUE),
      smc_none_q75 = quantile(smc_none_irr, 0.75, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # iii <- inci_wide %>% filter(as.Date(inci_wide$yearmonth) > as.Date('2017-06-01') &
  #                               as.Date(inci_wide$yearmonth) < as.Date('2018-01-01'))
  # Filter for specific date range
  iii_summary <- inci_summary %>% 
    filter(as.Date(yearmonth) > as.Date('2017-06-01') &
             as.Date(yearmonth) < as.Date('2018-01-01'))
  
  
  # Plot 1: Adding RTSS (filtered)
  ggplot(iii_summary) + 
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = both_smc_q25, ymax = both_smc_q75), 
                fill = '#FFE900', alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = both_smc_median, color = 'Both vs SMC'), linewidth = 1) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = rtss_none_q25, ymax = rtss_none_q75), 
                fill = '#E53D00', alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = rtss_none_median, color = 'RTSS vs none'), linewidth = 1) +
    ylim(c(0, 1)) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_date(breaks = '1 month',
                 labels = scales::label_date_short()) +
    scale_color_manual(values = c('Both vs SMC' = '#FFE900', 'RTSS vs none' = '#E53D00')) +
    labs(x = 'Date',
         y = '1-IRR',
         color = 'Comparison') + 
    theme_bw(base_size = 14)
  ggsave(paste0(path, outputsfolder,'/irr_addingrtss_filtered.pdf'), plot = last_plot())
  
  # Plot 2: Adding RTSS (non-filtered)
  ggplot(inci_summary) + 
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = both_smc_q25, ymax = both_smc_q75), 
                fill = '#FFE900', alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = both_smc_median, color = 'Both vs SMC'), linewidth = 1) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = rtss_none_q25, ymax = rtss_none_q75), 
                fill = '#E53D00', alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = rtss_none_median, color = 'RTSS vs none'), linewidth = 1) +
    ylim(c(0, 1)) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_date(breaks = '3 months',
                 labels = scales::label_date_short()) +
    scale_color_manual(values = c('Both vs SMC' = '#FFE900', 'RTSS vs none' = '#E53D00')) +
    labs(x = 'Date',
         y = '1-IRR',
         color = 'Comparison') + 
    theme_bw(base_size = 14)
  ggsave(paste0(path, outputsfolder,'/irr_addingrtss_nonfiltered.pdf'), plot = last_plot(), width = 16, height = 6)
  
  # Filter for SMC comparison
  iii_summary_smc <- inci_summary %>% 
    filter(as.Date(yearmonth) > as.Date('2017-06-01') & as.Date(yearmonth) < as.Date('2018-01-01'))
  
  # Plot 3: Adding SMC (filtered)
  ggplot(iii_summary_smc) + 
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = both_rtss_q25, ymax = both_rtss_q75), 
                fill = '#2ACBCB', alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = both_rtss_median, color = 'Both vs RTSS'), linewidth = 1) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = smc_none_q25, ymax = smc_none_q75), 
                fill = '#046865', alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = smc_none_median, color = 'SMC vs none'), linewidth = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_date(breaks = '1 month',
                 labels = scales::label_date_short()) +
    scale_color_manual(values = c('Both vs RTSS' = '#2ACBCB', 'SMC vs none' = '#046865')) +
    ylim(c(-0.5, 1)) +
    labs(x = 'Date',
         y = '1-IRR',
         color = 'Comparison') + 
    theme_bw(base_size = 14)
  ggsave(paste0(path, outputsfolder,'/irr_addingsmc_filtered.pdf'), plot = last_plot())
  
  # Plot 4: Adding SMC (non-filtered)
  ggplot(inci_summary) + 
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = both_rtss_q25, ymax = both_rtss_q75), 
                fill = '#2ACBCB', alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = both_rtss_median, color = 'Both vs RTSS'), linewidth = 1) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = smc_none_q25, ymax = smc_none_q75), 
                fill = '#046865', alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = smc_none_median, color = 'SMC vs none'), linewidth = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_date(breaks = '3 months',
                 labels = scales::label_date_short()) +
    scale_color_manual(values = c('Both vs RTSS' = '#2ACBCB', 'SMC vs none' = '#046865')) +
    ylim(c(-1, 1)) + 
    labs(x = 'Date',
         y = '1-IRR',
         color = 'Comparison') + 
    theme_bw(base_size = 14)
  ggsave(paste0(path, outputsfolder,'/irr_addingsmc_nonfiltered.pdf'), plot = last_plot(), width = 16, height = 6)
  
  # Print summary statistics
  cat("Mean of both_smc_irr (filtered):", mean(iii_summary$both_smc_median), "\n")
  cat("Mean of rtss_none_irr (filtered):", mean(iii_summary$rtss_none_median), "\n")
  cat("Mean of both_smc_irr (non-filtered):", mean(inci_summary$both_smc_median), "\n")
  cat("Mean of rtss_none_irr (non-filtered):", mean(inci_summary$rtss_none_median), "\n")
  
  cat("Mean of both_rtss_irr (filtered, Jun 2017-Jan 2018):", mean(iii_summary_smc$both_rtss_median), "\n")
  cat("Mean of smc_none_irr (filtered, Jun 2017-Jan 2018:", mean(iii_summary_smc$smc_none_median), "\n")
  cat("Mean of both_rtss_irr (non-filtered):", mean(inci_summary$both_rtss_median), "\n")
  cat("Mean of smc_none_irr (non-filtered):", mean(inci_summary$smc_none_median), "\n")
  
  # compare both vs smc to rtss vs none -- if synergy then could expect that the curve of both vs smc is shifted right 
  # ggplot(iii)+#nci_wide %>% filter(year <= 2018)) + 
  #   geom_line(aes(x = as.Date(yearmonth), y = both_smc_irr, group = sim_id, color = 'Both vs SMC')) +
  #   geom_line(aes(x = as.Date(yearmonth), y = rtss_none_irr, group = sim_id, color = 'RTSS vs none')) +
  #   geom_point(aes(x = as.Date(yearmonth), y = both_smc_irr, group = sim_id, color = 'Both vs SMC')) +
  #   geom_point(aes(x = as.Date(yearmonth), y = rtss_none_irr, group = sim_id, color = 'RTSS vs none')) +
  #   ylim(c(-0, 1)) +
  #   geom_hline(yintercept = 0, linetype = 2) +
  #   scale_x_date(breaks = '1 month',
  #                labels = scales::label_date_short()) +
  #   scale_color_manual(values = c('#FFE900', '#E53D00')) +
  #   labs(x = 'Date',
  #        y = '1-IRR',
  #        color = 'Comparison',
  #        caption = 'Each line is a simulation') + 
  #   theme_bw(base_size = 12)
  # ggsave(paste0(path, outputsfolder,'/irr_addingrtss_filtered.pdf'), plot = last_plot())
  # mean(iii$both_smc_irr)
  # mean(iii$rtss_none_irr)
  # 
  # ggplot(inci_wide) + 
  #   geom_line(aes(x = as.Date(yearmonth), y = both_smc_irr, group = sim_id, color = 'Both vs SMC')) +
  #   geom_line(aes(x = as.Date(yearmonth), y = rtss_none_irr, group = sim_id, color = 'RTSS vs none')) +
  #   geom_point(aes(x = as.Date(yearmonth), y = both_smc_irr, group = sim_id, color = 'Both vs SMC')) +
  #   geom_point(aes(x = as.Date(yearmonth), y = rtss_none_irr, group = sim_id, color = 'RTSS vs none')) +
  #   ylim(c(0, 1)) +
  #   geom_hline(yintercept = 0, linetype = 2) +
  #   scale_x_date(breaks = '1 month',
  #                labels = scales::label_date_short()) +
  #   scale_color_manual(values = c('#FFE900', '#E53D00')) +
  #   labs(x = 'Date',
  #        y = '1-IRR',
  #        color = 'Comparison',
  #        caption = 'Each line is a simulation') + 
  #   theme_bw(base_size = 12)
  # ggsave(paste0(path, outputsfolder,'/irr_addingrtss_nonfiltered.pdf'), plot = last_plot(), width = 16, height = 6)
  # mean(iii$both_smc_irr)
  # mean(iii$rtss_none_irr)
  # 
  # # compare both vs rtss to smc vs none -- could expect curve of both vs rtss is shifted right or higher for longer 
  # ggplot(iii %>% filter(as.Date(yearmonth) > as.Date('2017-07-01'))) + 
  #   geom_line(aes(x = as.Date(yearmonth), y = both_rtss_irr, group = sim_id, color = 'Both vs RTSS')) +
  #   geom_line(aes(x = as.Date(yearmonth), y = smc_none_irr, group = sim_id, color = 'SMC vs none')) +
  #   geom_point(aes(x = as.Date(yearmonth), y = both_rtss_irr, group = sim_id, color = 'Both vs RTSS')) +
  #   geom_point(aes(x = as.Date(yearmonth), y = smc_none_irr, group = sim_id, color = 'SMC vs none')) +
  #   geom_hline(yintercept = 0, linetype = 2) +
  #   scale_x_date(breaks = '1 month',
  #                labels = scales::label_date_short()) +
  #   scale_color_manual(values = c('#2ACBCB', '#046865')) +
  #   ylim(c(-0.2, 1)) + 
  #   labs(x = 'Date',
  #        y = '1-IRR',
  #        color = 'Comparison') + 
  #   theme_bw(base_size = 12)
  # ggsave(paste0(path, outputsfolder,'/irr_addingsmc_filtered.pdf'), plot = last_plot())
  # mean(iii[as.Date(iii$yearmonth) > as.Date('2017-07-01'),]$both_rtss_irr)
  # mean(iii[as.Date(iii$yearmonth) > as.Date('2017-07-01'),]$smc_none_irr)
  # 
  # # not filtering by date
  # ggplot(inci_wide) + 
  #   geom_line(aes(x = as.Date(yearmonth), y = both_rtss_irr, group = sim_id, color = 'Both vs RTSS')) +
  #   geom_line(aes(x = as.Date(yearmonth), y = smc_none_irr, group = sim_id, color = 'SMC vs none')) +
  #   geom_point(aes(x = as.Date(yearmonth), y = both_rtss_irr, group = sim_id, color = 'Both vs RTSS')) +
  #   geom_point(aes(x = as.Date(yearmonth), y = smc_none_irr, group = sim_id, color = 'SMC vs none')) +
  #   geom_hline(yintercept = 0, linetype = 2) +
  #   scale_x_date(breaks = '1 month',
  #                labels = scales::label_date_short()) +
  #   scale_color_manual(values = c('#2ACBCB', '#046865')) +
  #   ylim(c(-1, 1)) + 
  #   labs(x = 'Date',
  #        y = '1-IRR',
  #        color = 'Comparison') + 
  #   theme_bw(base_size = 12)
  # ggsave(paste0(path, outputsfolder,'/irr_addingsmc_nonfiltered.pdf'), plot = last_plot(), width = 16, height = 6)
  
}