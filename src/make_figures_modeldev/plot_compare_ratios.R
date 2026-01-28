# Script to compare the ratios of multiple timings 

plot_compare_ratios <- function(output_folders, cohort_folder = 'sim_cohort_generic'){
  
  path <- paste0('R:/Kelly/synergy_orderly/src/', cohort_folder, '/outputs/')
  
  # output_folders <- c('outputs_2026-01-19_3', # 122, 68
  #                     'outputs_2026-01-19_4', # 135, 80
  #                     'outputs_2026-01-19_2') # 100, 60
                      
  
  incidence_df <- purrr::map_dfr(
    output_folders,
    ~ readRDS(paste0(path, .x, "/incidence_summary.rds")) |>
      mutate(output_folder = .x)
  )
  
  incidence_df <- incidence_df %>%
    mutate(
      rtsstiming = case_when(
        output_folder == 'outputs_2026-01-19_3' ~ 68,
        output_folder == 'outputs_2026-01-19_4' ~ 80,
        output_folder == 'outputs_2026-01-19_2' ~ 60,
        output_folder == 'outputs_2026-01-23_15' ~ 68,
        output_folder == 'outputs_2026-01-23_18' ~ 55,
        output_folder == 'outputs_2026-01-23_19' ~ 55,
        output_folder == 'outputs_2026-01-23_20' ~ 80,
        output_folder == 'outputs_2026-01-23_21' ~ 80,
        TRUE ~ NA),
      smctiming = case_when(
        output_folder == 'outputs_2026-01-19_3' ~ 122,
        output_folder == 'outputs_2026-01-19_4' ~ 135,
        output_folder == 'outputs_2026-01-19_2' ~ 100,
        output_folder == 'outputs_2026-01-23_15' ~ 122,
        output_folder == 'outputs_2026-01-23_18' ~ 140,
        output_folder == 'outputs_2026-01-23_19' ~ 100,
        output_folder == 'outputs_2026-01-23_20' ~ 100,
        output_folder == 'outputs_2026-01-23_21' ~ 140,
        TRUE ~ NA)
    ) %>%
    mutate(scenario = case_when(
      output_folder == 'outputs_2026-01-23_15' ~ 'Balanced',
      output_folder == 'outputs_2026-01-23_19' ~ 'Early',
      output_folder == 'outputs_2026-01-23_21' ~ 'Late',
      output_folder == 'outputs_2026-01-23_18' ~ 'Early vaccine, late SMC',
      output_folder == 'outputs_2026-01-23_20' ~ 'Late vaccine, early SMC',
      TRUE ~ NA
    ))
  
  # metadata_df <- readRDS(paste0(path, outputsfolder, "/metadata_df.rds"))
  # base_inputs <- readRDS(paste0(path, outputsfolder, "/base_inputs.rds"))
  # params <- readRDS(paste0(path, outputsfolder, "/parameter_grid.rds"))
  # # smc_dates <- as.Date(unlist(formatted$smc_dose_days[10][1:4]), origin = '2017-04-01')
  # smc_dates <- readRDS('R:/Kelly/synergy_orderly/shared/median_smc_dates.rds') %>%
  #   filter(country == base_inputs$country) %>%
  #   pull(date)
  # # smc_dates <- as.Date(unlist(all$smc_dose_days[11][1:4]), origin = '2017-04-01')
  # smc_lines <- data.frame(
  #   xintercept = rep(smc_dates,2),
  #   arm = rep(c('smc', 'both'), each = length(smc_dates)),
  #   color = '#709176'
  # )
  # # metadata_df$vaccination_day[1] = 90
  # rtss_lines <- data.frame(
  #   xintercept = as.Date(rep(c(mean(metadata_df$vaccination_day)-60, mean(metadata_df$vaccination_day)-30, mean(metadata_df$vaccination_day), 
  #                              mean(metadata_df$vaccination_day)[1]+364, mean(metadata_df$vaccination_day)+730),2), origin = '2017-04-01'),
  #   arm = rep(c('rtss','both'), length(6)),
  #   color = '#59114D'
  # )
  
  
  # Plot of ratios of expected vs predicted efficacy of both vs none 
  ratioplot <- ggplot(incidence_df) + 
    # geom_vline(data = smc_lines, aes(xintercept = xintercept, color = 'SMC delivery'), 
    #            linetype = 2, linewidth = 0.8, alpha = 0.7)+
    # geom_vline(data = rtss_lines, aes(xintercept = xintercept, color = 'RTS,S delivery'), 
    #            linetype = 3, linewidth = 0.8, alpha = 0.7) +
    geom_line(aes(x = as.Date(yearmonth), y = ratio_pred_exp_median, group = scenario,#paste0('RTSS: ', rtsstiming, '; SMC: ', smctiming), 
                  color = scenario),#paste0('RTSS: ', rtsstiming, '; SMC: ', smctiming)), 
              # color = '#3E6990', 
              alpha = 0.7,
              linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = 2) +
    scale_x_date(breaks = '3 months',
                 labels = scales::label_date_short()) +
    scale_color_brewer(palette = 'Dark2') +
    # scale_color_manual(values = c('SMC delivery' = '#709176',
    #                               'RTS,S delivery' = '#470024')) +#'#449DD1'
    labs(x = 'Date',
         y = 'Ratio of modelled/expected',
         color = 'Scenario',
         fill = 'Scenario') + 
    theme_bw(base_size = 14)
  
  ratioplot
  ggsave('R:/Kelly/synergy_orderly/figures/predicted_vs_expected_compare_ratios.pdf', plot = ratioplot, width = 10, height = 4)
  
}

# plot_compare_ratios(output_folders = c('outputs_2026-01-19_3', # 122, 68
#                       'outputs_2026-01-19_4', # 135, 80
#                       'outputs_2026-01-19_2') # 100, 60
#                     )
  