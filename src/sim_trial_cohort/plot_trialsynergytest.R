# Function to plot the best-fitting parameters with the repetitions of full coverage/median timings 

plot_synergytest <- function(){
  cohort_folder <- 'sim_trial_cohort'
  path <- paste0('R:/Kelly/synergy_orderly/src/', cohort_folder, '/outputs/')
  # synergy test -- 2-arm fitting (before finished)
  # outputs_folders <-  c('BF' = 'outputs_2026-03-26_4',
  #                       'Mali' = 'outputs_2026-03-26_3'
  #                       )
   # 3-arm fitting -- synergy test 
  outputs_folders <-  c('BF' = 'outputs_2026-03-25_6',
                        'Mali' = 'outputs_2026-03-25_5'#
  )
  # # Folders below are from the bestreps runs for the 2-arm fitting (before finished)
  # outputs_folders <- c('Mali' = 'outputs_2026-03-26_2', 
  #                   'BF' = 'outputs_2026-03-26' ## 
  # )
  # # Folders below are from the bestreps runs for the 3-arm fitting 
  # outputs_folders <- c('Mali' = 'outputs_2026-03-24_2', # 
  #                   'BF' = 'outputs_2026-03-24' # -- 
  # )
  
  smc_lines_mali <- readRDS('R:/Kelly/synergy_orderly/shared/median_smc_dates.rds') %>%
    ungroup() %>%
    filter(country == 'Mali' & arm != 'rtss') %>%
    select(date, arm, country) %>%
    mutate(color = '#709176')
  smc_lines_bf <- readRDS('R:/Kelly/synergy_orderly/shared/median_smc_dates.rds') %>%
    ungroup() %>%
    filter(country == 'BF' & arm != 'rtss') %>%
    select(date, arm, country) %>%
    mutate(color = '#709176')
  smc_lines <- bind_rows(smc_lines_bf, smc_lines_mali)
  
  rtss_lines_mali <-  readRDS('R:/Kelly/synergy_orderly/shared/median_rtss_dates.rds') %>% ungroup() %>%
    filter(country == 'Mali' & arm != 'smc') %>%
    select(date, arm, country) %>% 
    mutate(color = '#59114D') 
  rtss_lines_bf <-  readRDS('R:/Kelly/synergy_orderly/shared/median_rtss_dates.rds') %>% ungroup() %>%
    filter(country == 'BF' & arm != 'smc') %>%
    select(date, arm, country) %>% 
    mutate(color = '#59114D') 
  rtss_lines <- bind_rows(rtss_lines_mali, rtss_lines_bf)
  colors <- RColorBrewer::brewer.pal(9, 'BuPu')
  colors <- c(colors[4:length(colors)], 'black')
  
  # pull in data for halfyear 
  inci_summary_halfyear_BF <- readRDS(paste0(path, outputs_folders['BF'], '/inci_summary_all_halfyear.rds')) %>%
    mutate(country = 'BF')
  inci_summary_halfyear_Mali <- readRDS(paste0(path, outputs_folders['Mali'], '/inci_summary_all_halfyear.rds')) %>%
    mutate(country = 'Mali')
  
  inci_summary_halfyear <- bind_rows(inci_summary_halfyear_BF, inci_summary_halfyear_Mali)
  
  # Pull in data for yearmonth
  inci_summary_yearmonth_BF <- readRDS(paste0(path, outputs_folders['BF'], '/inci_summary_wide_yearmonth.rds')) %>%
    mutate(country = 'BF')
  inci_summary_yearmonth_Mali <- readRDS(paste0(path, outputs_folders['Mali'], '/inci_summary_wide_yearmonth.rds')) %>%
    mutate(country = 'Mali')
  
  inci_summary_yearmonth <- bind_rows(inci_summary_yearmonth_BF, inci_summary_yearmonth_Mali) %>%
    rename(yearmonth = time_value)
  
  
  # Plot of the ratio of model-predicted to expected by aggregation unit 
  ggplot(inci_summary_halfyear %>% filter(metric == 'difference_inci_averted_pred_exp'), 
         aes(x = time_value, y = median, 
             color = time_value, group = time_value)) +
    # model estimated
    geom_point(size = 1) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                  width = 0.2, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = colors) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
    labs(
      x = NULL,#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Half-year',
      y = "Difference in model-predicted versus expected\ncases averted per 1000 people of\ncombination vs no intervention",
      shape = NULL, linetype = NULL,
      color = NULL#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Study half-year'
    ) +
    facet_wrap(~country) + 
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = 'none')
  
  ggsave(paste0('R:/Kelly/synergy_orderly/src/', cohort_folder, '/thesis_plots/difference_inci_average_byhalfyear_BFMali.pdf'), 
         plot = last_plot(), width = 9, height = 6)
  
  
  # Plot of expected vs predicted cases averted per 1000 of both vs none 
  exppred_inci <- ggplot(inci_summary_yearmonth %>% filter(yearmonth > '2017-05-01')) + 
    geom_vline(data = smc_lines, aes(xintercept = date, color = 'SMC delivery', group = country), linetype = 2, linewidth = 0.8)+
    geom_vline(data = rtss_lines, aes(xintercept = date, color = 'RTS,S delivery', group = country), linetype = 3, linewidth = 0.8) +
    geom_line(aes(x = as.Date(yearmonth), y = inci_averted_expected_median, color = 'Expected cases averted\nper 1000 person months'), linewidth = 1) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = inci_averted_expected_q025, ymax = inci_averted_expected_q975,
                    fill = 'Expected cases averted\nper 1000 person months'), alpha = 0.3) +
    geom_line(aes(x = as.Date(yearmonth), y = inci_averted_model_median, color = 'Model-predicted cases averted\nper 1000 person months'), 
              linewidth = 1) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = inci_averted_model_q025, ymax = inci_averted_model_q975,
                    fill = 'Model-predicted cases averted\nper 1000 person months'), alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_x_date(breaks = '3 months',
                 labels = scales::label_date_short()) +
    scale_color_manual(values = c('Model-predicted cases averted\nper 1000 person months' = '#59C9A5',
                                  'SMC delivery' = '#709176',
                                  'RTS,S delivery' = '#470024',
                                  'Expected cases averted\nper 1000 person months' = '#6457A6')) +#'#449DD1'
    scale_fill_manual(values = c('Expected cases averted\nper 1000 person months' = '#6457A6',
                                 'Model-predicted cases averted\nper 1000 person months' = '#59C9A5')) +#'#449DD1'
    facet_wrap(~country,
               nrow = 2) + 
    labs(x = 'Date',
         y = 'Cases averted per 1000 people',
         color = NULL,
         fill = NULL) + 
    theme_bw(base_size = 14)
  ggsave(paste0('R:/Kelly/synergy_orderly/src/', cohort_folder, '/thesis_plots/predicted_vs_expected_combined_incidenceaverted_BFMali.pdf'), 
         plot = exppred_inci, width = 12, height = 9)
  
  # # Plot of difference of expected vs predicted cases averted per 1000 of both vs none 
  # diffplot <-  ggplot(inci_summary_yearmonth %>% filter(yearmonth > '2017-05-01')) + 
  #   geom_vline(data = smc_lines, aes(xintercept = date, color = 'SMC delivery'), 
  #              linetype = 2, linewidth = 0.8, alpha = 0.7)+
  #   geom_vline(data = rtss_lines, aes(xintercept = date, color = 'RTS,S delivery'), 
  #              linetype = 3, linewidth = 0.8, alpha = 0.7) +
  #   geom_line(aes(x = as.Date(yearmonth), y = difference_inci_averted_pred_exp_median), 
  #             color = '#3E6990', linewidth = 1) +
  #   geom_ribbon(aes(x = as.Date(yearmonth), ymin = difference_inci_averted_pred_exp_q025, 
  #                   ymax = difference_inci_averted_pred_exp_q975),
  #               fill = '#3E6990', alpha = 0.5) +
  #   geom_hline(yintercept = 1, linetype = 2) +
  #   scale_x_date(breaks = '3 months',
  #                labels = scales::label_date_short()) +
  #   scale_color_manual(values = c('SMC delivery' = '#709176',
  #                                 'RTS,S delivery' = '#470024')) +#'#449DD1'
  #   # scale_fill_manual(values = c('Expected efficacy' = '#6457A6',
  #   #                              'Model-predicted efficacy' = '#59C9A5')) +#'#449DD1'
  #   labs(x = 'Date',
  #        y = 'Difference in model-predicted cases averted\nper 1000 people versus expected cases\naverted per 1000 people',
  #        color = NULL,
  #        fill = NULL) + 
  #   theme_bw(base_size = 14)
  # 
  # combineddifference <- plot_grid(exppred_inci, diffplot + theme(legend.position = 'none'), nrow = 2,
  #                                 align = 'v')
  # ggsave(paste0(path, outputsfolder,'/predicted_vs_expected_combined_incidence.pdf'), plot = diffplot, width = 10, height = 4)
  # ggsave(paste0(path, outputsfolder,'/predicted_vs_expected_combined_incidence_and_difference.pdf'), plot = combineddifference, width = 12, height = 9)
  
  
}