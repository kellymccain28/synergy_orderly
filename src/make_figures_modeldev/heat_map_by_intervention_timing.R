# FUnction to make heat map that compares the overall value of 'synergy' as defined by the impact difference in terms of 
# model-predicted - expected incidence averted 

#' folders is a vector of output folders in src/sim_cohort_generic/outputs/ which correspond to the timings we want to compare
#' extralabel is a character vector to modify file name from paste0(path, '/heat_plot_comparison', Sys.Date(), extralabel,'.pdf')
heat_map_by_intervention_timing <- function(folders = c('outputs_2026-02-18_9',
                                                        'outputs_2026-02-18_8',
                                                        'outputs_2026-02-18_7',
                                                        'outputs_2026-02-18_6',
                                                        'outputs_2026-02-18_5',
                                                        'outputs_2026-02-18_4',
                                                        'outputs_2026-02-18_3', # this is scenario 5 used in thesis ch 5
                                                        'outputs_2026-02-18_25',
                                                        'outputs_2026-02-18_24',
                                                        'outputs_2026-02-18_23',
                                                        'outputs_2026-02-18_22',
                                                        'outputs_2026-02-18_21',
                                                        'outputs_2026-02-18_20',
                                                        'outputs_2026-02-18_2',
                                                        'outputs_2026-02-18_19',
                                                        'outputs_2026-02-18_18',
                                                        'outputs_2026-02-18_17',
                                                        'outputs_2026-02-18_16',
                                                        'outputs_2026-02-18_15',
                                                        'outputs_2026-02-18_14',
                                                        'outputs_2026-02-18_13',
                                                        'outputs_2026-02-18_12',
                                                        'outputs_2026-02-18_11',
                                                        'outputs_2026-02-18_10',
                                                        'outputs_2026-02-18'),
                                            extralabel = NULL){
  cohort_folder = 'sim_cohort_generic'
  path <- paste0('R:/Kelly/synergy_orderly/src/', cohort_folder, '/outputs/')
  
  # Pull in data from all output folders and get the summarized value of the difference metric 
  overall_inci <- purrr::map_dfr(
    folders,
    ~ readRDS(paste0(path, .x, "/inci_summary_overall.rds")) |>
      mutate(output_folder = .x) %>%
      select(output_folder, time_value, contains('difference_cases_averted_pred_exp'), 
             contains('difference_inci_averted_pred_exp'))
  )
  
  smc_day <- purrr::map_dfr(
    folders,
    ~ readRDS(paste0(path, .x, "/parameter_grid.rds")) |>
      mutate(output_folder = .x) %>% ungroup() %>%
      mutate(first_smc_dose = smc_dose_days[[1]][1]) %>%
      slice_head() %>% select(output_folder, first_smc_dose)
  )
  
  vax_day <- purrr::map_dfr(
    folders,
    ~ readRDS(paste0(path, .x, "/metadata_df.rds")) |>
      mutate(output_folder = .x) %>% ungroup() %>%
      slice_head() %>% select(output_folder, vaccination_day)
  )
  
  overall_inci_summ <- overall_inci %>%
    left_join(smc_day) %>% 
    left_join(vax_day) %>% ungroup() %>%
    mutate(vaccination_day = as.Date(vaccination_day, origin = "2017-04-01"),
           first_smc_dose = as.Date(first_smc_dose, origin = "2017-04-01")) %>%
    mutate(# flag if the difference is statistically sig or not 
           significant = ifelse(difference_cases_averted_pred_exp_q025 > 0 | difference_cases_averted_pred_exp_q975 < 0, '*', ''))
  
  # Make the heat plot 
  ggplot(overall_inci_summ) +
    geom_tile(aes(x = vaccination_day, y = first_smc_dose, 
                  fill = difference_inci_averted_pred_exp_median)) + 
    geom_text(aes(x = vaccination_day, y = first_smc_dose,
                  label = significant)) +
    scale_fill_gradient2(low = '#F07167',
                         mid = 'white',
                         high = '#0081A7') + 
    labs(fill = 'Median difference in cases averted\nper 1000 person-months\n(predicted - expected)',
         x = 'Date of 3rd vaccine dose',
         y = 'Date of first monthly SMC dose') + 
    theme_classic(base_size = 14)
  
  ggsave(filename = paste0(path, 'heat_plot_comparison', Sys.Date(), extralabel,'.pdf'), 
         plot = last_plot(),
         height = 6, width = 12)
  
  # Calculate percentage difference
  overall_averted <- purrr::map_dfr(
    folders,
    ~ readRDS(paste0(path, .x, "/inci_summary_all_halfyear.rds")) |>
      mutate(output_folder = .x) %>%
      filter(metric == 'inci_averted_model' | metric == 'inci_averted_expected') %>%
      filter(time_value == 'Overall')
      
      # select(output_folder, time_value, contains('difference_cases_averted_pred_exp'), 
      #        contains('difference_inci_averted_pred_exp'))
  )
  overall_avertedjoin <- overall_averted %>%
    left_join(smc_day) %>% 
    left_join(vax_day) %>% ungroup() %>%
    mutate(vaccination_day = as.Date(vaccination_day, origin = "2017-04-01"),
           first_smc_dose = as.Date(first_smc_dose, origin = "2017-04-01"))
  
  inci_comparison <- overall_avertedjoin %>%
    select(time_value, time_value_num, metric, median, output_folder, vaccination_day, first_smc_dose) %>%
    pivot_wider(id_cols = c(time_value, time_value_num,output_folder, vaccination_day, first_smc_dose), 
                names_from = metric, 
                values_from = median) %>%
    mutate(
      pct_diff = (inci_averted_model - inci_averted_expected) / inci_averted_expected * 100,
      pct_diff_abs = abs(pct_diff),  # optional: absolute percentage difference
      diff_direction = ifelse(pct_diff > 0, "model higher", "model lower")
    )
  
  # Get heat plot with percent differences 
  ggplot(inci_comparison) +
    geom_tile(aes(x = vaccination_day, y = first_smc_dose, 
                  fill = pct_diff)) + 
    # geom_text(aes(x = vaccination_day, y = first_smc_dose, 
    #               label = significant)) +
    scale_fill_gradient2(low = '#F07167',
                         mid = 'white',
                         high = '#0081A7') + 
    labs(fill = 'Median percent difference in cases averted\nper 1000 person-months\n(predicted - expected)',
         x = 'Date of 3rd vaccine dose',
         y = 'Date of first monthly SMC dose') + 
    theme_classic(base_size = 14)
  
  ggsave(filename = paste0(path, 'heat_plot_comparison_perc', Sys.Date(), extralabel,'.pdf'), 
         plot = last_plot(),
         height = 6, width = 12)
  
  saveRDS(inci_comparison, paste0(path, '/inci_percent_averted.rds'))
}
# 
# heat_map_by_intervention_timing(folders = c('outputs_2026-03-23_5',
#                                             'outputs_2026-03-23_6',
#                                             'outputs_2026-03-23_7',
#                                             'outputs_2026-03-23_8',
#                                             'outputs_2026-03-23_9',
#                                             'outputs_2026-03-23_10',
#                                             'outputs_2026-03-23_11',
#                                             'outputs_2026-03-23_12',
#                                             'outputs_2026-03-23_13'),
#                                 extralabel = '_time_immunity')
# 
# heat_map_by_intervention_timing(folders = c('outputs_2026-03-24_9',
#                                             'outputs_2026-03-24_8',
#                                             'outputs_2026-03-24_7',
#                                             'outputs_2026-03-24_6',
#                                             'outputs_2026-03-24_5',
#                                             'outputs_2026-03-24_4',
#                                             'outputs_2026-03-24_3',
#                                             'outputs_2026-03-24_2',
#                                             'outputs_2026-03-24'
# ),
#                                 extralabel = '_time_immunity_100')
# 
# heat_map_by_intervention_timing(folders = c('outputs_2026-03-24_18',
#                                             'outputs_2026-03-24_17',
#                                             'outputs_2026-03-24_16',
#                                             'outputs_2026-03-24_15',
#                                             'outputs_2026-03-24_14',
#                                             'outputs_2026-03-24_13',
#                                             'outputs_2026-03-24_12',
#                                             'outputs_2026-03-24_11',
#                                             'outputs_2026-03-24_10'),
# extralabel = '_time_immunity_100_ANDsm')
