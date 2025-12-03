# Plot efficacy by year
plot_efficacy <- function(outputsfolder){
  # Load packages
  library(zoo)
  library(survival)
  library(survminer)
  library(broom)
  library(ggplot2)
  library(tidyverse)
  
  source("R:/Kelly/synergy_orderly/shared/get_cox_efficacy.R")
  source("R:/Kelly/synergy_orderly/shared/likelihood.R")
  
  path <- 'R:/Kelly/synergy_orderly/src/sim_cohort_generic/outputs/'
  # outputsfolder <- 'outputs_2025-12-01_2'
  
  # Using the outputs from monthly_incidence_plot.R
  formatted <- readRDS(paste0(path, outputsfolder, '/formatted_infrecords.rds'))
  inci <- readRDS(paste0(path, outputsfolder, '/incidence.rds'))
  
  # Get efficacy by year and overall
  smcrefresults <- get_cox_efficacy(df = formatted,
                                    ref = 'arm_smcref',
                                    model = TRUE)
  
  
  ## RTSS comparator by year and overall ----
  rtssrefresults <- get_cox_efficacy(df = formatted,
                                     ref = 'arm_rtssref',
                                     model = TRUE)
  
  #None comparator by year and overall ----
  nonerefresults <- get_cox_efficacy(df = formatted,
                                     ref = 'arm_noneref',
                                     model = TRUE)
  
  # Plot the vaccine efficacies
  tidy_results <- rbind(smcrefresults,
                        rtssrefresults,
                        nonerefresults) %>%
    filter(term !='None vs. SMC' & term != 'None vs. RTSS') %>%
    mutate(reference = case_when(
      grepl('vs. SMC', term) ~ 'SMC reference',
      grepl('vs. RTSS', term) ~ 'RTS,S reference',
      grepl('vs None', term) ~ 'No intervention\nreference',
      TRUE ~ NA
    ),
    term = factor(term, levels = c("RTSS vs None", "SMC vs None", "Both vs None",
                                   "SMC vs. RTSS", "Both vs. RTSS",
                                   "RTSS vs. SMC", "Both vs. SMC")),
    year = factor(year, 
                  levels = c(1, 2, 3, 'overall'),
                  labels = c("Year 1", "Year 2", "Year 3", "Overall")))
  # filter(term %in% c("Both vs. RTSS",'RTSS vs. SMC','Both vs. SMC'))
  
  # expected:
  rtssnone <- tidy_results[tidy_results$term == 'RTSS vs None', c('VE', "VE_lower", "VE_upper")]
  smcnone <- tidy_results[tidy_results$term == 'SMC vs None',c('VE', "VE_lower", "VE_upper")]#$VE
  years <- tidy_results[tidy_results$term == 'RTSS vs None',]$year
  
  expected_efficacy <- rtssnone * (1 - smcnone) + smcnone
  expected_efficacy$year <- tidy_results[tidy_results$term == 'RTSS vs None',]$year
  
  tidy_results <- tidy_results %>%
    left_join(expected_efficacy %>%
                rename(VE_expected = VE, 
                       VE_lower_expected = VE_lower, 
                       VE_upper_expected = VE_upper), by = 'year') %>%
    mutate(VE_expected = ifelse(term == 'Both vs None', VE_expected, NA),
           VE_lower_expected = ifelse(term == 'Both vs None', VE_lower_expected, NA),
           VE_upper_expected = ifelse(term == 'Both vs None', VE_upper_expected, NA))
  # tidy_results %>% filter(year == 'Overall', reference == 'No intervention\nreference') 
    
  # 0.44892964*(1-0.53744008) +0.53744008
  # [1] 0.7450969  -- synergy this time! 12-01_2 with updated SMC -- what changed is the smc parameters (#2), population size 
  # # this should be lower than the both vs none if there is synergy
  
  efficacy_plot <- ggplot(tidy_results %>% filter(year == 'Overall')) +
    geom_point(aes(x = term, y = VE, color = reference), size = 1.5) +
    geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper, color = reference), width = 0.1, linewidth = 0.8) +
    geom_point(aes(x = term, y = VE_expected, color = 'Expected efficacy'), 
               size = 1.5) +
    geom_errorbar(aes(x = term, ymin = VE_lower_expected, ymax = VE_upper_expected, color = 'Expected efficacy'), 
                  width = 0.1, linewidth = 0.8) +
    geom_hline(aes(yintercept = 0)) +
    labs(y = "Efficacy",
         x = NULL, 
         color = 'Comparison group') +
    scale_color_manual(values = c(
      'SMC reference' = '#2292A4',
      'RTS,S reference' = '#BDBF09',
      'No intervention\nreference' = '#D96C06',
      'Expected efficacy' = colorspace::lighten('#D96C06', amount = 0.4))) + 
    scale_y_continuous(breaks = seq(min(floor(tidy_results$VE_lower * 10)/10), 1.1, 0.2),
                       limits = c(min(ceiling(tidy_results$VE_lower * 10)/10),1),
                       labels = scales::label_percent()) +
    # facet_wrap(~ year) + 
    theme_bw(base_size = 12) 
  ggsave(paste0(path, outputsfolder,'/model_efficacy.pdf'), plot = last_plot(), height = 6, width = 10)
  # efficacy_plot
  
  return(tidy_results)
}
