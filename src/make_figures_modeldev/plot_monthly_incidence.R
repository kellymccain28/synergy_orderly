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
  
  inci <- bind_rows(inci_all)
  saveRDS(inci, paste0(path, outputsfolder, '/incidence.rds'))
  
  files_formatted <- files[grepl('formatted', files)]
  files_formatted <- files_formatted[grepl('batch', files_formatted)]
  formatted_all <- lapply(files_formatted, readRDS)
  
  all <- bind_rows(formatted_all)
  rm(formatted_all)
  saveRDS(all, paste0(path, outputsfolder, '/formatted_infrecords.rds'))
  
  metadata_df <- readRDS(paste0(path, outputsfolder, "/metadata_df.rds"))
  base_inputs <- readRDS(paste0(path, outputsfolder, "/base_inputs.rds"))
  params <- readRDS(paste0(path, outputsfolder, "/parameter_grid.rds"))
  
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
    if(!is.null(unlist(all$smc_dose_days[1]))){
      smc_dates <- as.Date(unlist(all$smc_dose_days[1]), origin = '2017-04-01')
    } else if(!is.null(unlist(all$smc_dose_days[2]))){
      smc_dates <- as.Date(unlist(all$smc_dose_days[2]), origin = '2017-04-01')
    } else if(!is.null(unlist(all$smc_dose_days[3]))){
      smc_dates <- as.Date(unlist(all$smc_dose_days[3]), origin = '2017-04-01')
    } else if(!is.null(unlist(all$smc_dose_days[4]))){
      smc_dates <- as.Date(unlist(all$smc_dose_days[4]), origin = '2017-04-01')
    } else if(!is.null(unlist(all$smc_dose_days[5]))){
      smc_dates <- as.Date(unlist(all$smc_dose_days[5]), origin = '2017-04-01')
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
  ggplot(inci_summ)+
    # geom_vline(xintercept = as.Date(unlist(all$smc_dose_days[1][1:4]), origin = '2017-04-01'), linetype = 2, color = '#4D9DE0') +
    geom_vline(data = smc_lines, aes(xintercept = xintercept, color = 'SMC delivery'), linetype = 2, linewidth = 0.7)+
    geom_vline(data = rtss_lines, aes(xintercept = xintercept, color = 'RTS,S delivery'), linetype = 3, linewidth = 0.8) +
    geom_line(aes(x = as.Date(yearmonth), y = incidence_per_1000pm, color = arm), linewidth = 0.8) +
    # geom_point(aes(x = as.Date(yearmonth), y = incidence_per_1000pm/100, color = arm), size = 4) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = incidence_per_1000pm_lower, ymax = incidence_per_1000pm_upper, fill = arm),# color = arm),
                  linewidth = 1, width = 25, alpha = 0.3) +
    # geom_line(data = pbite, aes(x = date, y = prob_lagged*100), alpha = 0.2) +
    scale_x_date(breaks = '3 months',
                 labels = scales::label_date_short()) +
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
  
}
