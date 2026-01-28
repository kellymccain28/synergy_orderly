# Script to produce nice outputs of the monthly incidence of the generic cohort 
plot_monthly_incidence <- function(outputsfolder, cohort_folder = 'sim_cohort_generic'){
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
  # outputsfolder <- 'outputs_2025-12-01_2'
  
  # sim_results <- readRDS(paste0(path, outputsfolder, "/sim_results.rds"))
  files <- list.files(paste0(path, outputsfolder),
                      full.names = TRUE)
  files_inci <- files[grepl('incidence', files)]
  files_inci <- files_inci[grepl('batch', files_inci)]
  inci_all <- lapply(files_inci, readRDS)
  
  files_formatted <- files[grepl('formatted', files)]
  files_formatted <- files_formatted[grepl('batch', files_formatted)]
  formatted_all <- lapply(files_formatted, readRDS)
  
  # infectionrecords <-  readRDS(paste0(path, outputsfolder, "/infection_records.rds")) %>%
  #   # remove any infections that occurred within 7 days 
  #   group_by(rid, sim_id) %>%
  #   arrange(rid, sim_id, detection_day) %>%
  #   mutate(previous_detday = lag(detection_day),
  #          diff = detection_day - previous_detday) %>%
  #   filter(diff > 7 | is.na(diff)) %>% select(-diff, -previous_detday)
  metadata_df <- readRDS(paste0(path, outputsfolder, "/metadata_df.rds"))
  base_inputs <- readRDS(paste0(path, outputsfolder, "/base_inputs.rds"))
  params <- readRDS(paste0(path, outputsfolder, "/parameter_grid.rds"))
  
  # infectionrecords <- purrr::map_df(sim_results, "infection_records")
  # params <- purrr::map_df(sim_results, 'params')
  
  # Format for incidence calculation 
  # formatted_all <- lapply(params$sim_id, function(x){
  #   format_model_output(model_data = infectionrecords,
  #                       cohort = 'generic',
  #                       simulation = x)})
  all <- bind_rows(formatted_all)
  saveRDS(all, paste0(path, outputsfolder, '/formatted_infrecords.rds'))
  
  # formattedinfrecords_monthly <- lapply(params$sim_id, function(x){
  #   format_model_output_monthly(model_data = infectionrecords, 
  #                               cohort = 'generic',
  #                               simulation = x, 
  #                               n_months = 36)
  # })
  # allmonths <- bind_rows(formattedinfrecords_monthly)
  
  # Calculate incidence 
  # inci_all <- lapply(params$sim_id, function(x){
  #   aa <- all %>% filter(sim_id == x)
  #   
  #   get_incidence(df_children = metadata_df,
  #                 casedata = aa) %>%
  #     mutate(sim_id = x)
  # })
  inci <- bind_rows(inci_all)
  saveRDS(inci, paste0(path, outputsfolder, '/incidence.rds'))
  
  
  # Summarize over all simulations 
  inci_summ <- inci %>%
    group_by(arm, year, month, yearmonth) %>%
    summarize(across(c(incidence_per_1000pm),
                     list(lower = ~quantile(.x, 0.025, na.rm = TRUE),
                          median = ~quantile(.x, 0.5, na.rm = TRUE),
                          upper = ~quantile(.x, 0.975, na.rm = TRUE)),
                     .names = "{.col}_{.fn}") ) %>%
    # rename those variables with _median to be just the variable name 
    rename_with(.fn = \(x)sub("_median","", x)) 
  
  # Plot incidence with p_bite overlaid 
  pbite <- params$p_bite[[1]]
  if(cohort_folder == 'sim_cohort_generic'){
    smc_dates <- as.Date(unlist(formatted$smc_dose_days[11][1:4]), origin = '2017-04-01')
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
  ggplot(inci_summ)+
    # geom_vline(xintercept = as.Date(unlist(all$smc_dose_days[1][1:4]), origin = '2017-04-01'), linetype = 2, color = '#4D9DE0') +
    geom_vline(data = smc_lines, aes(xintercept = xintercept, color = 'SMC delivery'), linetype = 2, linewidth = 0.7)+
    geom_vline(data = rtss_lines, aes(xintercept = xintercept, color = 'RTS,S delivery'), linetype = 3, linewidth = 0.8) +
    geom_line(aes(x = as.Date(yearmonth), y = incidence_per_1000pm, color = arm), linewidth = 0.8) +
    # geom_point(aes(x = as.Date(yearmonth), y = incidence_per_1000pm/100, color = arm), size = 4) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = incidence_per_1000pm_lower, ymax = incidence_per_1000pm_upper, fill = arm),# color = arm),
                  linewidth = 1, width = 25, alpha = 0.3) +
    # geom_line(data = pbite, aes(x = date, y = prob_lagged*100), alpha = 0.2) +
    scale_color_manual(values =  c('both' = '#E15554', 
                                   'none' = '#E1BC29',
                                   'rtss' = '#3BB273',
                                   'smc' = '#7768AE',
                                   'SMC delivery' = '#709176',
                                   'RTS,S delivery' = '#470024'),
                       breaks = c('both','none','rtss','smc','SMC delivery','RTS,S delivery'))+#c('#C44536','#772E25','#197278','#283D3B'))+
    scale_fill_manual(values =  c('both' = '#E15554', 
                                   'none' = '#E1BC29',
                                   'rtss' = '#3BB273',
                                   'smc' = '#7768AE',
                                   'SMC delivery' = '#709176',
                                   'RTS,S delivery' = '#470024'),
                       breaks = c('both','none','rtss','smc','SMC delivery','RTS,S delivery'))+#c('#C44536','#772E25','#197278','#283D3B'))+
    labs(color = 'Intervention arm',
         x = 'Date',
         y = 'Incidence per 1000 person-months') +
    theme_bw(base_size = 14) + 
    facet_wrap(~factor(arm, levels = c('none','rtss','smc','both')))
  
  ggsave(paste0(path, outputsfolder,'/model_monthly_incidence.pdf'), plot = last_plot(),
         height = 8, width = 14)
  
  
  
  # Plot efficacy over time using incidence curves 
  # none <- inci %>%
  #   filter(arm == 'none') %>%
  #   rename(incidence_per_1000pm_none = incidence_per_1000pm,
  #          rate_none = rate, 
  #          person_months_none = person_months) %>% 
  #   select(sim_id, year, month, yearmonth, incidence_per_1000pm_none, rate_none, person_months_none)
  # rtss <- inci %>%
  #   filter(arm == 'rtss') %>%
  #   rename(incidence_per_1000pm_rtss = incidence_per_1000pm,
  #          rate_rtss = rate, 
  #          person_months_rtss = person_months) %>% 
  #   select(sim_id, year, month, yearmonth, incidence_per_1000pm_rtss, rate_rtss, person_months_rtss)
  # smc <- inci %>%
  #   filter(arm == 'smc') %>%
  #   rename(incidence_per_1000pm_smc = incidence_per_1000pm,
  #          rate_smc = rate, 
  #          person_months_smc = person_months) %>% 
  #   select(sim_id, year, month, yearmonth, incidence_per_1000pm_smc, rate_smc, person_months_smc)
  # 
  # df <- inci %>%
  #   filter(arm == 'both') %>%
  #   left_join(none) %>%
  #   left_join(rtss) %>%
  #   left_join(smc) %>%
  #   # Calculate efficacy
  #   mutate(rtss_none_eff = 1 - (incidence_per_1000pm_rtss / incidence_per_1000pm_none),
  #          smc_none_eff = 1 - (incidence_per_1000pm_smc / incidence_per_1000pm_none),
  #          both_none_eff = 1 - (incidence_per_1000pm / incidence_per_1000pm_none),
  #          both_rtss_eff = 1 - (incidence_per_1000pm / incidence_per_1000pm_rtss),
  #          both_smc_eff = 1 - (incidence_per_1000pm / incidence_per_1000pm_smc)) %>%
  #   pivot_longer(cols = rtss_none_eff:both_smc_eff,
  #                names_to = 'efficacy_type',
  #                values_to = 'efficacy')
  # 
  # ggplot(df %>% filter(sim_id == 'parameter_set_6_generic_0.9')) + 
  #   geom_line(aes(x = yearmonth, y = efficacy, color = efficacy_type)) + ylim(c(0,1))
}
