# Function to plot the best-fitting parameters with the repetitions of full coverage/median timings 
# Function to plot the syntest and bestreps comparison for the 3-arm fit 
# this is for chapter 6 (synergy in fitted model section) to show that coverage is
# not driving the synergy metric 
# adapted from the 1-IRR average script
plot_synergytest <- function(){
  cohort_folder <- 'sim_trial_cohort'
  path <- paste0('R:/Kelly/synergy_orderly/src/', cohort_folder, '/outputs/')
  # synergy test -- 2-arm fitting (before finished)
  # outputs_folders <-  c('BF' = 'outputs_2026-03-26_4',
  #                       'Mali' = 'outputs_2026-03-26_3'
  #                       )
  # 3-arm fitting -- synergy test 
  outputs_folders <-  c('BF3syn' = 'outputs_2026-03-25_6',
                        'Mali3syn' = 'outputs_2026-03-25_5',
                        'BF2syn' = 'outputs_2026-03-30_7',
                        'Mali2syn' = 'outputs_2026-03-30_8',
                       
                        'BF3best' = 'outputs_2026-03-24',
                        'Mali3best' = 'outputs_2026-03-24_2',
                        'BF2best' = 'outputs_2026-03-30',
                        'Mali2best' = 'outputs_2026-03-30_2')
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
    dplyr::select(date, arm, country) %>%
    mutate(color = '#709176')
  smc_lines_bf <- readRDS('R:/Kelly/synergy_orderly/shared/median_smc_dates.rds') %>%
    ungroup() %>%
    filter(country == 'BF' & arm != 'rtss') %>%
    dplyr::select(date, arm, country) %>%
    mutate(color = '#709176')
  smc_lines <- bind_rows(smc_lines_bf, smc_lines_mali)
  
  rtss_lines_mali <-  readRDS('R:/Kelly/synergy_orderly/shared/median_rtss_dates.rds') %>% ungroup() %>%
    filter(country == 'Mali' & arm != 'smc') %>%
    dplyr::select(date, arm, country) %>% 
    mutate(color = '#59114D') 
  rtss_lines_bf <-  readRDS('R:/Kelly/synergy_orderly/shared/median_rtss_dates.rds') %>% ungroup() %>%
    filter(country == 'BF' & arm != 'smc') %>%
    dplyr::select(date, arm, country) %>% 
    mutate(color = '#59114D') 
  rtss_lines <- bind_rows(rtss_lines_mali, rtss_lines_bf)
  colors <- c('syntest' = '#3ABB94',
              'bestreps' = '#FB9D4B')
  
  # pull in data for halfyear ( these were saved from summarize_IRRs.R )
  inci_summary_halfyear_BFsyn3 <- readRDS(paste0(path, outputs_folders['BF3syn'], '/inci_summary_all_halfyear.rds')) %>%
    mutate(country = 'BF') %>%
    mutate(type = 'syntest')
  inci_summary_halfyear_Malisyn3 <- readRDS(paste0(path, outputs_folders['Mali3syn'], '/inci_summary_all_halfyear.rds')) %>%
    mutate(country = 'Mali') %>%
    mutate(type = 'syntest')
  inci_summary_halfyear_BFbest3 <- readRDS(paste0(path, outputs_folders['BF3best'], '/inci_summary_all_halfyear.rds')) %>%
    mutate(country = 'BF') %>%
    mutate(type = 'bestreps')
  inci_summary_halfyear_Malibest3 <- readRDS(paste0(path, outputs_folders['Mali3best'], '/inci_summary_all_halfyear.rds')) %>%
    mutate(country = 'Mali') %>%
    mutate(type = 'bestreps')
  
  inci_summary_halfyear3 <- bind_rows(inci_summary_halfyear_BFsyn3, inci_summary_halfyear_Malisyn3,
                                      inci_summary_halfyear_BFbest3, inci_summary_halfyear_Malibest3) %>%
    mutate(armfit = '3-arm fit')
  
  inci_summary_halfyear_BFsyn2 <- readRDS(paste0(path, outputs_folders['BF2syn'], '/inci_summary_all_halfyear.rds')) %>%
    mutate(country = 'BF')%>%
    mutate(type = 'syntest')
  inci_summary_halfyear_Malisyn2 <- readRDS(paste0(path, outputs_folders['Mali2syn'], '/inci_summary_all_halfyear.rds')) %>%
    mutate(country = 'Mali')%>%
    mutate(type = 'syntest')
  inci_summary_halfyear_BFbest2 <- readRDS(paste0(path, outputs_folders['BF2best'], '/inci_summary_all_halfyear.rds')) %>%
    mutate(country = 'BF') %>%
    mutate(type = 'bestreps')
  inci_summary_halfyear_Malibest2 <- readRDS(paste0(path, outputs_folders['Mali2best'], '/inci_summary_all_halfyear.rds')) %>%
    mutate(country = 'Mali') %>%
    mutate(type = 'bestreps')
  
  inci_summary_halfyear2 <- bind_rows(inci_summary_halfyear_BFsyn2, inci_summary_halfyear_Malisyn2,
                                      inci_summary_halfyear_BFbest2, inci_summary_halfyear_Malibest2) %>%
    mutate(armfit = '2-arm fit')
  
  inci_summary_halfyear <- bind_rows(inci_summary_halfyear2, inci_summary_halfyear3) %>%
    filter(metric == 'difference_inci_averted_pred_exp')
  
  # Plot of the difference of model-predicted to expected between bestreps for 2 and 3-arm fits of BF and Mali
  ggplot(inci_summary_halfyear, 
         aes(x = time_value, y = median, 
             color = type, group = type)) +
    # model estimated
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
    geom_point(size = 1, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                  width = 0.5, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = colors,
                       labels = c('syntest' = 'Full coverage',
                                  'bestreps' = 'Trial coverage')) +
    # Add annotations for synergistic/antagonistic
    annotate("text", x = -Inf, y = max(inci_summary_halfyear$upper_ci[inci_summary_halfyear$metric == 'difference_inci_averted_pred_exp'], na.rm = TRUE) * 0.03, 
             label = "Synergistic", color = "green4", hjust = -0.1, vjust = 0, size = 3, fontface = "italic") +
    annotate("text", x = -Inf, y = min(inci_summary_halfyear$lower_ci[inci_summary_halfyear$metric == 'difference_inci_averted_pred_exp'], na.rm = TRUE) * 0.03, 
             label = "Antagonistic", color = "red3", hjust = -0.1, vjust = 1, size = 3, fontface = "italic") +
    labs(
      x = NULL,#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Half-year',
      y = "Difference in model-predicted versus expected\ncases averted per 1000 person-months of\ncombination vs no intervention",
      shape = NULL, linetype = NULL,
      color = NULL#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Study half-year'
    ) +
    facet_grid(row = vars(country), col = vars(armfit),
               scales = 'free') + 
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(paste0('R:/Kelly/synergy_orderly/src/', cohort_folder, '/thesis_plots/difference_inci_average_byhalfyear_BFMali_2and3.pdf'), 
         plot = last_plot(), width = 12, height = 10)
  
  # Plot of the difference of model-predicted to expected between bestreps for 3-arm fit of BF and Mali
  ggplot(inci_summary_halfyear %>% filter(armfit == '3-arm fit'), 
         aes(x = time_value, y = median, 
             color = type, group = type)) +
    # model estimated
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
    geom_point(size = 1, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                  width = 0.5, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = colors,
                       labels = c('syntest' = 'Full coverage',
                                  'bestreps' = 'Trial coverage')) +
    # Add annotations for synergistic/antagonistic
    annotate("text", x = -Inf, y = max(inci_summary_halfyear$upper_ci[inci_summary_halfyear$metric == 'difference_inci_averted_pred_exp'], na.rm = TRUE) * 0.03, 
             label = "Synergistic", color = "green4", hjust = -0.1, vjust = 0, size = 3, fontface = "italic") +
    annotate("text", x = -Inf, y = min(inci_summary_halfyear$lower_ci[inci_summary_halfyear$metric == 'difference_inci_averted_pred_exp'], na.rm = TRUE) * 0.03, 
             label = "Antagonistic", color = "red3", hjust = -0.1, vjust = 1, size = 3, fontface = "italic") +
    labs(
      x = NULL,#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Half-year',
      y = "Difference in model-predicted versus expected\ncases averted per 1000 person-months of\ncombination vs no intervention",
      shape = NULL, linetype = NULL,
      color = NULL#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Study half-year'
    ) +
    facet_wrap(~country) + 
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(paste0('R:/Kelly/synergy_orderly/src/', cohort_folder, '/thesis_plots/difference_inci_average_byhalfyear_BFMali_3armfit.pdf'), 
         plot = last_plot(), width = 11, height = 6)
  
  return(inci_summary_halfyear %>% filter(time_value == 'Overall' & armfit == '3-arm fit'))
  
}
plot_synergytest()
