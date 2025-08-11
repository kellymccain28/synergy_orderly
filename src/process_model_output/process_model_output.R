# Replicate survival analysis/calculation of incidence by month of trial results to model outputs
library(grates) # ifw e want to use greates_isoweek format like incidence does
library(survival)
library(tidyverse)
library(broom)
library(survminer)

orderly_strict_mode()

orderly_dependency(name = "sim_cohort_grid", "latest()",
                   files = c('simulation_outputs/'))
                   # files = c("outputs/infection_records.rds", 
                   #           "outputs/parasitemia.rds", 
                   #           "outputs/metadata_children.rds",
                   #           "inputs.rds"))

orderly_resource(files = c('format_model_output.R',
                           'read_in_outputs.R'))
source('format_model_output.R')
source('read_in_outputs.R')

orderly_shared_resource('get_cox_efficacy.R')
source('get_cox_efficacy.R')

orderly_shared_resource('get_incidence.R')
source('get_incidence.R')

dir.create('outputs/')
dir.create("outputs/plots/")

# orderly_artefact(description = 'plots and formatted model and analysis output',
#                  files = c(
#                    'outputs/plots/efficacy_plot.png',
#                    'outputs/model_output_formatted.rds',
#                    'outputs/plots/monthly_incidence_model.png',
#                    'outputs/surv_analysis_model.rds', 
#                    'outputs/monthly_incidence_model.rds'
#                  ))

# infection_records <- readRDS("outputs/infection_records.rds")
# parasitemia <- readRDS("outputs/parasitemia.rds")
# metadata_child <- readRDS("outputs/metadata_children.rds")
# inputs <- readRDS("inputs.rds") # this is just to carry over the inputs 

# Read in results from the latest run of sim_cohort_grid
results <- read_in_outputs()

baseline_inputs <- results[[1]]
infection_records <- results[[2]]$infection_records
parasitemia <- results[[2]]$parasitemia_data
parameters <- results[[2]]$parameters
metadata_child <- results[[2]]$child_metadata

sim_ids <- parameters$sim_id

# Format model output 
outputs <- lapply(sim_ids, function(sim){
  format_model_output(infection_records, simulation = sim)
})
saveRDS(outputs, 'outputs/model_outputs_formatted.rds')



analyse_output <- function(outputs, simulation){
  
  # Make the directories to save outputs
  sim_dir <- file.path('outputs', simulation)
  
  if (!dir.exists(sim_dir)) {
    dir.create(sim_dir, recursive = TRUE)
  }
  
  plots_dir <- file.path(sim_dir, "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir)
  }
  
  # filter to individual sim ids
  output <- outputs %>%
    bind_rows() %>%
    filter(sim_id == simulation)
  
  # Get number of cases and person time per month 
  monthly <- output %>%
    mutate(month = lubridate::month(detection_day)) %>%
    group_by(arm, month)
  
  # Do the same analysis as we do on trial data in trial_results.R
  model_output <- output %>%
    mutate(arm_smcref = factor(arm, levels = c('smc','vax','vaxsmc','none')),
           arm_rtssref = factor(arm, levels = c('vax','smc','vaxsmc','none')),
           arm_noneref = factor(arm, levels = c('none','vax','smc','vaxsmc')),
           country = 'BF') %>%
    filter(start_time != end_time)
  
  
  # Survival analysis to reproduce results from trial
  ## SMC comparator by year and overall ---- 
  smcrefresults <- get_cox_efficacy(df = model_output, 
                                    ref = 'arm_smcref',
                                    model = TRUE)
  
  ## RTSS comparator by year and overall ---- 
  rtssrefresults <- get_cox_efficacy(df = model_output, 
                                     ref = 'arm_rtssref',
                                     model = TRUE)
  
  #None comparator by year and overall ---- 
  nonerefresults <- get_cox_efficacy(df = model_output, 
                                     ref = 'arm_noneref',
                                     model = TRUE)
  
  # Plot the vaccine efficacies 
  tidy_results <- rbind(smcrefresults, 
                        rtssrefresults, 
                        nonerefresults)
  
  saveRDS(tidy_results, paste0('outputs/', simulation, '/surv_analysis_model.rds'))
  
  
  # Plot results from survival analysis 
  survresults <- tidy_results %>%
    filter(term %in% c("Both vs. RTSS",'RTSS vs. SMC','Both vs. SMC'))
  
  efficacy_plot <- ggplot(survresults)+
    geom_point(aes(x = term, y = VE)) + 
    geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper), width = 0.2) +
    geom_hline(aes(yintercept = 0)) + 
    labs(y = "Vaccine efficacy",
         x = NULL) +
    scale_y_continuous(breaks = seq(min(floor(survresults$VE_lower * 10)/10), 1, 0.1),
                       limits = c(min(floor(survresults$VE_lower * 10)/10),1)) + 
    facet_wrap(~ year)
  
  ggsave(filename = paste0('outputs/', simulation, '/plots/efficacy_plot.png'), efficacy_plot, height = 6, width = 6)
  
  # Get expected efficacy of combination for each year ----
  # Calculate expected vs observed
  ve_comparison <- tidy_results %>%
    group_by(year) %>%
    filter(term %in% c("SMC vs None", "RTSS vs None", "Both vs None")) %>%
    select(term, VE) %>%
    pivot_wider(names_from = term, values_from = VE) %>%
    mutate(
      expected_both = 1 - (1 - `SMC vs None`) * (1 - `RTSS vs None`),
      observed_both = `Both vs None`,
      difference = observed_both - expected_both,
      interaction_type = case_when(
        difference > 0 ~ "Synergistic",
        difference < 0 ~ "Antagonistic", 
        TRUE ~ "Independent"
      )
    )
  
  saveRDS(ve_comparison, file = paste0('outputs/', simulation, '/expected_efficacies.rds'))
  
  # Get monhtly incidence
  monthly_inci_model <- get_incidence(model = TRUE, 
                                      df_children = metadata_child, 
                                      casedata = model_output)
  
  saveRDS(monthly_inci_model, paste0('outputs/', simulation, '/monthly_incidence_model.rds'))
  # for fitting, can use date, person_months, n_cases intervention variables 
  
  incidence_plot <- monthly_inci_model  %>%
    # filter(intervention != 'none') %>%
    ggplot(aes(x = date, y = incidence_per_1000pm, color = arm)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower_per_1000, ymax = upper_per_1000, color = arm),
                  alpha = 0.9, width = 10, linewidth = 1) +
    geom_line(linewidth = 1) + 
    # geom_hline(aes(xintercept = c())) +
    scale_x_date(breaks = '3 months',
                 labels = scales::label_date_short()) + 
    # facet_wrap(~arm, nrow = 4) +
    # scale_y_continuous(breaks = seq(0,150,25)) +
    labs(
      title = "Monthly malaria incidence per 1000 person-months",
      x = "Year",
      y = "Incidence (per 1000 person-months)",
      color = "Study arm",
      fill = "Study arm"
    ) +
    theme_minimal(base_size = 16)
  
  
  ggsave(filename = paste0('outputs/', simulation, '/plots/monthly_incidence_model.png'), bg = 'white', incidence_plot, height = 6, width = 12)
  
  
  # Number of infections per day # incidence
  ints <- c('smc','vax','vaxsmc','none')
  weeks <- seq.Date(as.Date('2017-04-02'), as.Date('2020-03-31'), 7)
  allweeks <- expand.grid(arm = ints, 
                          week = weeks)# %>%
  # mutate(week = rep(1:length(weeks), each = length(ints)))
  
  weeklyinf <- infection_records %>%
    mutate(week = floor_date(as.Date(detection_day, origin = '2017-04-01'), unit = 'week')) %>%
    merge( allweeks, all.y = TRUE) %>%
    group_by(arm, week) %>%
    count()
  
  n_infections_weekly <- ggplot(weeklyinf) + 
    geom_line(aes(x = week, y = n, color = arm, group = arm)) +
    geom_point(aes(x = week, y = n, color = arm, group = arm)) +
    facet_wrap(~ arm) +
    labs(y = 'Weekly number of infections',
         x = 'Year since start of follow up period',
         caption = 'Assuming that the liver stage lasts 8 days') +
    theme(legend.position = 'none') + 
    theme_bw(base_size = 15)
  
  ggsave(filename = paste0('outputs/', simulation, '/plots/n_infections_weekly_model.png'), bg = 'white', 
         n_infections_weekly, height = 6, width = 12)
  
  n_infections_monthly <- ggplot(cases_by_month_model) + 
    geom_point(aes(x = monthyear, y = n_cases, color = arm, group = arm)) +
    geom_line(aes(x = monthyear, y = n_cases, color = arm, group = arm)) +
    facet_wrap(~ arm) +
    labs(y = 'Monthly number of infections',
         x = 'Year since start of follow up period',
         caption = 'Assuming that the liver stage lasts 8 days') +
    theme(legend.position = 'none') + 
    theme_bw(base_size = 15)
  
  ggsave(filename = paste0('outputs/', simulation, '/plots/n_infections_monthly_model.png'), bg = 'white', 
         n_infections_monthly, height = 6, width = 12)
  
  ggplot(parasitemia %>% filter(day1_BSinfection == 19 & arm == 'smc')) + 
    geom_line(aes(x=time_ext/7, y = prob_smckill))
  
  # Proportion detectable ---- 
  prop_det <- ggplot(infection_records ) + 
    geom_bar(aes(x = arm, group = as.factor(detectable), fill = as.factor(detectable)),
             position = 'fill') + 
    scale_fill_manual(values = c('darkmagenta','goldenrod')) +
    labs(#title = 'Detectable infections in each arm group',
      fill = 'Detectable')
  
  ggsave(filename = paste0('outputs/', simulation, '/plots/prop_detectable_model.png'), prop_det, height = 6, width = 6)
  
    # Cumulative incidence ---- 
  output <- output %>%
    mutate(end_time = ifelse(start_time==end_time, end_time + 0.0001, end_time),
           arm = factor(arm, levels = c('vaxsmc','vax','smc','none')))
  
  kmsurvobj <- survfit(Surv(start_time, 
                            end_time, 
                            event) ~ arm,#+ strata(country), 
                       data = output)
  
  cum_inci <- ggsurvplot(kmsurvobj, #group.by = "country",
                         # linetype = "strata",
                         # facet.by = 'country',
                         tables.theme = theme_cleantable(),
                         conf.int = TRUE,
                         fun = 'cumhaz',
                         # risk.table = TRUE,
                         # cumevents = TRUE,
                         ggtheme = theme_bw(base_size = 16),
                         palette = 'Dark2',
                         censor.size = 3)
  
  ggsave(filename = paste0('outputs/', simulation, '/plots/cum_inci_model.png'), cum_inci$plot, height = 6, width = 10)
  
  
  survival <- ggsurvplot(kmsurvobj, #group.by = "country",
                         # linetype = "strata",
                         # facet.by = 'country',
                         tables.theme = theme_cleantable(),
                         conf.int = TRUE,
                         # fun = 'cumhaz',
                         # risk.table = TRUE,
                         # cumevents = TRUE,
                         ggtheme = theme_bw(base_size = 16),
                         palette = 'Dark2',
                         censor.size = 3)
  
  ggsave(filename = paste0('outputs/', simulation, '/plots/survival_model.png'), survival$plot, height = 6, width = 10)
  
  message("Saved outputs for ", simulation)
}

lapply(sim_ids, function(sim){
  analyse_output(outputs, simulation = sim)
  }  )
