# Script to run all figure-producing scripts that analyse the outputs of the cohort simulations
setwd('R:/Kelly/synergy_orderly/src/make_figures_modeldev')

filestosource <- list.files(full.names = TRUE)
source("./plot_1-IRR.R")
source('./plot_1-IRR_average.R')
source("./plot_monthly_incidence.R")
source("./plot_initial_merozoites.R")
source('./plot_hazard_ratios.R')
source('./plot_time_to_threshold.R')
source('./plot_compare_ratios.R')

outputs_folder <- "outputs_2025-12-02_2_treat_1start_137threshold5000"
outputs_folder <- "outputs_2025-12-02_3_treat_0.9start_150threshold5000"
outputs_folder <- "outputs_2025-12-02_treat_0start_137threshold5000"
outputs_folder <- "outputs_2025-12-02_treat_0.9start_137threshold10000"
outputs_folder <- "outputs_2025-12-05_treat_0.9start_40threshold5000"
# outputs_folder <- "outputs_2025-12-08_treat_0.9start_141threshold5000" #used for init mero plot
outputs_folder <- "outputs_2025-12-08_treat_0.9start_140threshold5000"
outputs_folder <- "outputs_2025-12-08_treat_0.9start_122threshold5000" # **  folder used in thesis (before 20 dec), seasonal, but with ab+tcell pars
outputs_folder <- "outputs_2025-12-09_treat_0.9start_115threshold5000_perennial"
outputs_folder <- "outputs_2025-12-09_treat_0.9start_115threshold5000_noboosters"
outputs_folder <- "outputs_2025-12-11_treat_0.9start_115threshold5000_noboost_constant15"
outputs_folder <- "outputs_2025-12-11_treat_0.9start_115threshold5000_constant15_notdiv5"
outputs_folder <- 'outputs_2025-12-17_treat_0.9start_115threshold5000' # ** changing to this one; seasonal generic, ab only parameters
outputs_folder <- 'outputs_2025-12-17_treat_0.9start_115threshold5000_2' # perennial generic, ab only parameters
outputs_folder <- 'outputs_2025-12-17_treat_0.9start_115threshold5000_3' # perennial, ab only parameters, vacc day 75
outputs_folder <- 'outputs_2025-12-18_treat_0.9start_115threshold5000' # pperennial, ab only parameters, earlier vaccination (day 20)
outputs_folder <- 'outputs_2025-12-21_treat_0.9start_122threshold5000' # seasonal, ab only, 75 vaccination, 

#Updated RTSS pars (15 Jan 2026) 

# perennial:
outputs_folder <- 'outputs_2026-01-15'
outputs_folder <- 'outputs_2026-01-22'
# constant: 
outputs_folder <- 'outputs_2026-01-15_2'


# seasonal: 
outputs_folder <- 'outputs_2026-01-15_3'

outputs_folder <- 'outputs_2026-01-16_4'
outputs_folder <- 'outputs_2026-01-16_5'
outputs_folder <- 'outputs_2026-01-16_6'
outputs_folder <- 'outputs_2026-01-19'
outputs_folder <- 'outputs_2026-01-19_2'
outputs_folder <- 'outputs_2026-01-19_3' # *used in thesis as of 20 Jan
outputs_folder <- 'outputs_2026-01-19_4'
outputs_folder <- 'outputs_2026-01-19_6'
outputs_folder <- 'outputs_2026-01-19_7'

# perennial (1.66, 3.45, 0.00311) 
outputs_folder <- 'outputs_2026-01-23_8' #- done

# seasonal (1.66, 3.45, 0.00311) 
outputs_folder <- 'outputs_2026-01-23_9' # 100 80 - done
outputs_folder <- 'outputs_2026-01-23_10' # 100 80 - done
outputs_folder <- 'outputs_2026-01-23_11' # 100 80
outputs_folder <- 'outputs_2026-01-23_12' # 100 80
outputs_folder <- 'outputs_2026-01-23_13' # 100 80

# perennial (1.74, 4.69, 0.00259)
outputs_folder <- 'outputs_2026-01-23_16' # 122 68 - done

# seasonal (1.74, 4.69, 0.00259)
outputs_folder <- 'outputs_2026-01-23_15' # balanced 122 68 - done
outputs_folder <- 'outputs_2026-01-23_18' # early vax late smc 140 55 - done
outputs_folder <- 'outputs_2026-01-23_19' # early 100 55 - done
outputs_folder <- 'outputs_2026-01-23_20' # late vax early smc 100 80 - done
outputs_folder <- 'outputs_2026-01-23_21' # late 140 80

outputs_folder <- 'outputs_2026-01-26' # balanced with 32*3 runs

# Real trial simulations 
outputs_folder <- 'outputs_2026-01-27_8'


# first, do monthly inci, then order doesn't matter
plot_monthly_incidence(outputsfolder = outputs_folder)
plot_monthly_incidence(outputsfolder = outputs_folder,
                       cohort_folder = 'sim_trial_cohort')

# time to threshold
plot_time_to_threshold(outputsfolder = outputs_folder)
plot_time_to_threshold(outputsfolder = outputs_folder,
                       cohort_folder = 'sim_trial_cohort')

# efficacy - saved figure in folder and outputs efficacy results 
# hr_results <- plot_hazard_ratios(outputsfolder = outputs_folder)

# plot 1- IRR
plot_irr(outputsfolder = outputs_folder)
plot_irr(outputsfolder = outputs_folder,
         cohort_folder = 'sim_trial_cohort')

# plot 1-IRR average 
plot_irr_average(outputsfolder = outputs_folder, agg_unit = 'year')
plot_irr_average(outputsfolder = outputs_folder, agg_unit = 'halfyear')
plot_irr_average(outputsfolder = outputs_folder, agg_unit = 'year',
                 cohort_folder = 'sim_trial_cohort')
plot_irr_average(outputsfolder = outputs_folder, agg_unit = 'halfyear',
                 cohort_folder = 'sim_trial_cohort')


# compare ratios (only run after all folders in list have had monthly inci calculated, as this requires formatting )
plot_compare_ratios(output_folders = c('outputs_2026-01-23_15',
                                       'outputs_2026-01-23_18',
                                       'outputs_2026-01-23_19',
                                       'outputs_2026-01-23_20',
                                       'outputs_2026-01-23_21'))

# initial merozoites (only if we export parasitemia)
plot_initial_merozoites(outputsfolder = 'outputs_2026-01-23_22' ) #'outputs_2026-01-16', 'outputs_2025-12-08_treat_0.9start_141threshold5000', 'outputs_2026-01-15_4'
