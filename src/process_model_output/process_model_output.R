# Replicate survival analysis/calculation of incidence by month of trial results to model outputs
library(grates) # ifw e want to use greates_isoweek format like incidence does
library(survival)
library(tidyverse)
library(broom)
library(survminer)
library(orderly2)

orderly_strict_mode()

orderly_dependency(name = "sim_cohort_grid", "latest()",
                   files = c('simulation_outputs/'))
                   # files = c("outputs/infection_records.rds", 
                   #           "outputs/parasitemia.rds", 
                   #           "outputs/metadata_children.rds",
                   #           "inputs.rds"))

orderly_dependency(name = 'trial_results',
                   "latest()",
                   c('monthly_incidence_trial.rds'))

# key <- cyphr::data_key()
monthly_inci_trial <- readRDS("monthly_incidence_trial.rds")

orderly_resource(files = c('format_model_output.R',
                           'read_in_outputs.R',
                           'analyse_model_output.R',
                           'likelihood.R'))
source('format_model_output.R')
source('read_in_outputs.R')
source('likelihood.R')
source('analyse_model_output.R')

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

# baseline parameters same for all runs
baseline_inputs <- results[[1]] 
# linelist of people who were bitten by an infectious mosquito 
infection_records <- results[[2]]$infection_records
# parasitemia trajectories for each person 
# parasitemia <- results[[2]]$parasitemia_data
# parameters for each model run -- all are different
parameters <- results[[2]]$parameters
# information about interventions, etc for children in model 
metadata_child <- results[[2]]$child_metadata

sim_ids <- parameters$sim_id


# Format model output 
outputs <- lapply(sim_ids, function(sim){
  format_model_output(infection_records, simulation = sim)
})
saveRDS(outputs, 'outputs/model_outputs_formatted.rds')



# Initialise new column for the LL
parameters$ll <- NA

# parameters <- lapply(sim_ids, function(sim){
#   analyse_model_output(outputs, simulation = sim, parameters)
#   }  )
parameters <- Map(function(sim){
  analyse_model_output(outputs, simulation = sim, parameters)},
  sim_ids) 

parameters <- bind_rows(parameters)

saveRDS(parameters, 'outputs/parameters_ll.rds')

# plots comparing the different parameters and the log-likelihoods
# ggplot(parameters) +
#   geom_point(aes(x = max_SMC_kill_rate, y = ll, size = SMC_decay, color = season_start_day))
# ggplot(parameters) +
#   geom_point(aes(x = SMC_decay, y = ll, size = max_SMC_kill_rate, color = season_start_day))
llplot <- ggplot(parameters) + 
  geom_point(aes(x = season_start_day, y = ll, size = SMC_decay, color = max_SMC_kill_rate)) + 
  theme_bw()
ggsave('ll_plot.png', llplot)