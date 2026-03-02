# Script to compare the ratios of multiple timings 

plot_compare_ratios <- function(output_folders, cohort_folder = 'sim_cohort_generic'){
  
  path <- paste0('R:/Kelly/synergy_orderly/src/', cohort_folder, '/outputs/')
  
  incidence_df <- purrr::map_dfr(
    output_folders,
    ~ readRDS(paste0(path, .x, "/inci_summary_all_yearmonth.rds")) |>
      mutate(output_folder = .x) %>%
      rename(yearmonth = time_value) %>%
      filter(metric == 'difference_inci_averted_pred_exp')
  )
  
  smc_day <- purrr::map_dfr(
    output_folders,
    ~ readRDS(paste0(path, .x, "/parameter_grid.rds")) |>
      mutate(output_folder = .x) %>% ungroup() %>%
      mutate(first_smc_dose = smc_dose_days[[1]][1]) %>%
      slice_head() %>% select(output_folder, first_smc_dose)
  )
  
  vax_day <- purrr::map_dfr(
    output_folders,
    ~ readRDS(paste0(path, .x, "/metadata_df.rds")) |>
      mutate(output_folder = .x) %>% ungroup() %>%
      slice_head() %>% select(output_folder, vaccination_day)
  )
  
  incidence_df <- incidence_df %>%
    left_join(vax_day) %>%
    left_join(smc_day) %>%
    mutate(scenario = case_when(
      output_folder == 'outputs_2026-01-23_15' | output_folder == 'outputs_2026-01-26' ~ 'Balanced',
      output_folder == 'outputs_2026-01-23_19' ~ 'Early',
      output_folder == 'outputs_2026-01-23_21' ~ 'Late',
      output_folder == 'outputs_2026-01-23_18' ~ 'Early vaccine, late SMC',
      output_folder == 'outputs_2026-01-23_20' ~ 'Late vaccine, early SMC',
      
      output_folder == 'outputs_2026-02-10_2' ~ 'Balanced',
      output_folder == 'outputs_2026-02-10_5' ~ 'Early',
      output_folder == 'outputs_2026-02-10_6' ~ 'Late',
      output_folder == 'outputs_2026-02-10_8' ~ 'Early vaccine, late SMC',
      output_folder == 'outputs_2026-02-10_7' ~ 'Late vaccine, early SMC',
      TRUE ~ paste0('3rd vac dose: ', vaccination_day, ' 1st SMC dose: ', first_smc_dose)
    ),
    vaccination_day = paste0('3rd dose: ', as.Date(vaccination_day, origin = "2017-04-01")),
    first_smc_dose = as.Date(first_smc_dose, origin = "2017-04-01"))
  
  # Plot of ratios of expected vs predicted efficacy of both vs none 
  diffplot <- ggplot(incidence_df) + 
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = lower_ci, ymax = upper_ci,
                fill = as.factor(first_smc_dose)), alpha = 0.2) +
    geom_line(aes(x = as.Date(yearmonth), y = median, group = as.factor(first_smc_dose),
                  color = as.factor(first_smc_dose)),
              alpha = 0.8,
              linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = 2) +
    scale_x_date(breaks = '3 months',
                 labels = scales::label_date_short()) +
    scale_color_manual(values = colorRampPalette(c("#F4A261", "#2A9D8F"))(5),
                       labels = format(sort(unique(incidence_df$first_smc_dose)), "%d %b %Y")) +
    scale_fill_manual(values = colorRampPalette(c("#F4A261", "#2A9D8F"))(5),
                      labels = format(sort(unique(incidence_df$first_smc_dose)), "%d %b %Y")) +# scale_color_brewer(palette = 'Dark2') +
    # scale_fill_brewer(palette = 'Dark2') +
    # coord_cartesian(ylim = c(0.0,1.2))+
    labs(x = 'Date',
         y = 'Difference of cases averted per 1000 people\n(model-predicted - expected)',
         color = 'Date of first SMC dose',
         fill = 'Date of first SMC dose') + 
    theme_bw(base_size = 14) + 
    facet_wrap(~vaccination_day)
  
  diffplot
  ggsave(paste0('R:/Kelly/synergy_orderly/figures/predicted_vs_expected_compare_difference',Sys.Date(),'.pdf'),
         plot = diffplot, width = 14, height = 7)
  
}

# plot_compare_ratios(output_folders = c('outputs_2026-01-19_3', # 122, 68
#                       'outputs_2026-01-19_4', # 135, 80
#                       'outputs_2026-01-19_2') # 100, 60
#                     )
  