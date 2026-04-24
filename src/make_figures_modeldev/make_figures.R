# Script to run all figure-producing scripts that analyse the outputs of the cohort simulations
setwd('R:/Kelly/synergy_orderly/src/make_figures_modeldev')

library(purrr)
library(tidyr)
library(hipercow)
library(tidyverse)
# hipercow_environment_create(name = 'makeplots',
#                             sources = c("src/make_figures_modeldev/plot_1-IRR.R",
#                                         "src/make_figures_modeldev/plot_1-IRR_average.R",
#                                         "src/make_figures_modeldev/plot_monthly_incidence.R",
#                                         'src/make_figures_modeldev/plot_hazard_ratios.R',
#                                         'src/make_figures_modeldev/plot_time_to_threshold.R',
#                                         'src/make_figures_modeldev/plot_compare_ratios.R',
#                                         'src/make_figures_modeldev/bootstrap_metric.R',
#                                         'src/make_figures_modeldev/summarize_IRRs.R',
#                                         'src/make_figures_modeldev/make_figures_func.R'
#                             ))
# hipercow_provision(method = 'pkgdepends', environment = 'makeplots')

filestosource <- list.files(full.names = TRUE)
source("./plot_1-IRR.R")
source('./plot_1-IRR_average.R')
source("./plot_monthly_incidence.R")
source("./plot_initial_merozoites.R")
source('./plot_hazard_ratios.R')
source('./plot_time_to_threshold.R')
source('./plot_compare_ratios.R')
source('./bootstrap_metric.R')
source('./summarize_IRRs.R')


# Real trial simulations -- need to use cohort_folder = 'sim_trial_cohort'
cohort_folder = 'sim_trial_cohort'
outputs_folder <- 'outputs_2026-01-27_8'
outputs_folder <- 'outputs_2026-02-13_3'
outputs_folder <- 'outputs_2026-02-25_3'
outputs_folder <- 'outputs_2026-03-02' # redoing mali (no int arm had 0 inci before bc rid error)
outputs_folder <- 'outputs_2026-03-17_12'
outputs_folder <- 'outputs_2026-03-17_11'

# With updated SMC parameters fitted to threshold of 3000 -- as of 10 feb
cohort_folder <- 'sim_cohort_generic'
outputs_folders <- c('outputs_2026-02-10_2', # balanced 122 68 
                    'outputs_2026-02-10_4', # perennial 122 68
                    'outputs_2026-02-10_5', # early 100 55
                    'outputs_2026-02-10_6', # late 140 80
                    'outputs_2026-02-10_7', # late vax early smc 100 80
                    'outputs_2026-02-10_8') # early vax late smc 140 55

outputs_folders <- c('outputs_2026-02-12_8',
                    'outputs_2026-02-12_7',
                    'outputs_2026-02-12_6',
                    'outputs_2026-02-12_5',
                    'outputs_2026-02-12_4',
                    'outputs_2026-02-12_3',
                    'outputs_2026-02-12',
                    'outputs_2026-02-16',
                    'outputs_2026-02-17',
                    'outputs_2026-02-17_2',
                    'outputs_2026-02-17_3')

outputs_folders <- c('outputs_2026-02-17_6',
                     'outputs_2026-02-17_5',
                     'outputs_2026-02-17_4',
                     'outputs_2026-02-17_3',
                     'outputs_2026-02-17_2',
                     'outputs_2026-02-17')

# range of timings to use for heat map 
outputs_folders <- c('outputs_2026-02-18_9',
                     'outputs_2026-02-18_8',
                     'outputs_2026-02-18_7',
                     'outputs_2026-02-18_6',
                     'outputs_2026-02-18_5',
                     'outputs_2026-02-18_4',
                     'outputs_2026-02-18_3', # this is what is currently used in thesis for main results (as of 10 April, scenario 5)
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
                     'outputs_2026-02-18')


# no rtss decay
cohort_folder <- 'sim_cohort_generic'
outputs_folders <- c('outputs_2026-02-27',
                     paste0('outputs_2026-02-27_',seq(2, 9)))

# higher clearance threshold
cohort_folder <- 'sim_cohort_generic'
outputs_folders <- paste0('outputs_2026-02-27_',seq(9, 18))

# st immmunity, no gen adaptive immunity; 64 reps
cohort_folder <- 'sim_cohort_generic'
outputs_folders <- c('outputs_2026-03-20', #150,30
                     'outputs_2026-03-20_2'# 110,75
                     )

# generic cohrots to get parasitaemia over tiem; 64 reps
cohort_folder <- 'sim_cohort_generic'
outputs_folders <- c('outputs_2026-03-17', #110,75
                     'outputs_2026-03-17_2'# 150,30
)

outputs_folders <- c('outputs_2026-03-24_2','outputs_2026-03-24')
outputs_folders <- c('outputs_2026-03-23_5',
                    'outputs_2026-03-23_6',
                    'outputs_2026-03-23_7',
                    'outputs_2026-03-23_8',
                    'outputs_2026-03-23_9',
                    'outputs_2026-03-23_10',
                    'outputs_2026-03-23_11',
                    'outputs_2026-03-23_12',
                    'outputs_2026-03-23_13')

task_create_expr({
  cohort_folder <- 'sim_trial_cohort'
  outputs_folders <-  c('outputs_2026-03-17_4',
                        'outputs_2026-03-17_3',
                        'outputs_2026-03-24_2',
                        'outputs_2026-03-24',
                        'outputs_2026-03-30',
                        'outputs_2026-03-30_2')
  cohort_folder <- 'sim_cohort_generic'
  outputs_folders <- 'outputs_2026-02-18_3'
  outputs_folders <- c(#'outputs_2026-02-18_9',
                       # 'outputs_2026-02-18_8',
                       # 'outputs_2026-02-18_7',
                       # 'outputs_2026-02-18_6',
                       # 'outputs_2026-02-18_5',
                       # 'outputs_2026-02-18_4',
                       # # 'outputs_2026-02-18_3', # this is what is currently used in thesis for main results (as of 10 April, scenario 5)
                       # 'outputs_2026-02-18_25',
                       # 'outputs_2026-02-18_24',
                       # 'outputs_2026-02-18_23',
                       # 'outputs_2026-02-18_22',
                       # 'outputs_2026-02-18_21',
                       # 'outputs_2026-02-18_20',
                       # 'outputs_2026-02-18_2',
                       # 'outputs_2026-02-18_19',
                       # 'outputs_2026-02-18_18',
                       # 'outputs_2026-02-18_17',
                       # 'outputs_2026-02-18_16',
                       # 'outputs_2026-02-18_15',
                       # 'outputs_2026-02-18_14',
                       # 'outputs_2026-02-18_13',
                       # 'outputs_2026-02-18_12',
                       # 'outputs_2026-02-18_11',#*didn't finish, borke at plot_irr
                       'outputs_2026-02-18_10',
                       'outputs_2026-02-18')
  
  
  for(outputs_folder in outputs_folders){
    print(outputs_folder)
    # first, do monthly inci, then order doesn't matter
    plot_monthly_incidence(outputsfolder = outputs_folder,
                           cohort_folder = cohort_folder)
    
    message('monthly incidence plotted')
    
    # next, calculate the IRRs over simulations
    summarize_IRRs(outputsfolder = outputs_folder,
                   agg_unit = 'year')
    summarize_IRRs(outputsfolder = outputs_folder,
                   agg_unit = 'halfyear')
    summarize_IRRs(outputsfolder = outputs_folder,
                   agg_unit = 'yearmonth')
    message('finished summarising IRRs')
    
    # time to threshold
    plot_time_to_threshold(outputsfolder = outputs_folder,
                           cohort_folder = cohort_folder)
    message('plotted time to threshold')
    # efficacy - saved figure in folder and outputs efficacy results
    # hr_results <- plot_hazard_ratios(outputsfolder = outputs_folder)
    
    # plot 1- IRR
    plot_irr(outputsfolder = outputs_folder,
             cohort_folder = cohort_folder)
    message('plotted 1-IRR')
    
    # plot 1-IRR average
    plot_irr_average(outputsfolder = outputs_folder, agg_unit = 'year',
                     cohort_folder = cohort_folder)
    plot_irr_average(outputsfolder = outputs_folder, agg_unit = 'halfyear',
                     cohort_folder = cohort_folder)
    
    message('plotted 1-IRR averaged')
  }
  },
  environment = 'makeplots'
)

task_log_show('629dba2d5d4a24a54093db2d22e3d206') 

# compare ratios (only run after all folders in list have had monthly inci calculated, as this requires formatting )
plot_compare_ratios(output_folders = c('outputs_2026-02-18_9',
                                         'outputs_2026-02-18_8',
                                         'outputs_2026-02-18_7',
                                         'outputs_2026-02-18_6',
                                         'outputs_2026-02-18_5',
                                         'outputs_2026-02-18_4',
                                         'outputs_2026-02-18_3',
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
                                         'outputs_2026-02-18')) # late smc, early vax
# plot_compare_ratios(output_folders = c(#'outputs_2026-01-23_15',
#                                        'outputs_2026-01-23_18',
#                                        'outputs_2026-01-23_19',
#                                        'outputs_2026-01-23_20',
#                                        'outputs_2026-01-23_21',
#                                        'outputs_2026-01-26'))
# output_folders = c('outputs_2026-02-10_2', # balanced
#                    'outputs_2026-02-10_5', #early
#                    'outputs_2026-02-10_6', # late
#                    'outputs_2026-02-10_7', # early smc, late vax
#                    'outputs_2026-02-10_8')

heat_map_by_intervention_timing(folders = c('outputs_2026-02-18_9',
                                            'outputs_2026-02-18_8',
                                            'outputs_2026-02-18_7',
                                            'outputs_2026-02-18_6',
                                            'outputs_2026-02-18_5',
                                            'outputs_2026-02-18_4',
                                            'outputs_2026-02-18_3',
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
                                            'outputs_2026-02-18'),# all from 02-18
                                extralabel = '')

# initial merozoites (only if we export parasitemia)
plot_initial_merozoites(outputsfolder =  'outputs_2026-02-20_2') #'outputs_2026-01-23_22' 'outputs_2026-01-16', 'outputs_2025-12-08_treat_0.9start_141threshold5000', 'outputs_2026-01-15_4'
