# Script to run all figure-producing scripts that analyse the outputs of the cohort simulations
setwd('R:/Kelly/synergy_orderly/src/make_figures_modeldev')

filestosource <- list.files(full.names = TRUE)
source("./plot_1-IRR.R")
source("./plot_monthly_incidence.R")
source("./plot_initial_merozoites.R")
source('./plot_hazard_ratios.R')
source('./plot_time_to_threshold.R')

outputs_folder <- "outputs_2025-12-02_2_treat_1start_137threshold5000"
outputs_folder <- "outputs_2025-12-02_3_treat_0.9start_150threshold5000"
outputs_folder <- "outputs_2025-12-02_treat_0start_137threshold5000"
outputs_folder <- "outputs_2025-12-02_treat_0.9start_137threshold10000"
outputs_folder <- "outputs_2025-12-05_treat_0.9start_40threshold5000"
# outputs_folder <- "outputs_2025-12-08_treat_0.9start_141threshold5000"
outputs_folder <- "outputs_2025-12-08_treat_0.9start_140threshold5000"
outputs_folder <- "outputs_2025-12-08_treat_0.9start_122threshold5000" # ** current folder used in thesis 
outputs_folder <- "outputs_2025-12-09_treat_0.9start_115threshold5000_perennial"
outputs_folder <- "outputs_2025-12-09_treat_0.9start_115threshold5000_noboosters"
outputs_folder <- "outputs_2025-12-11_treat_0.9start_115threshold5000_noboost_constant15"
outputs_folder <- "outputs_2025-12-11_treat_0.9start_115threshold5000_constant15_notdiv5"
outputs_folder <- 'outputs_2025-12-17_treat_0.9start_115threshold5000' # seasonal generic, ab only parameters
outputs_folder <- 'outputs_2025-12-17_treat_0.9start_115threshold5000_2' # perennial generic, ab only parameters

# first, do monthly inci, then order doesn't matter
plot_monthly_incidence(outputsfolder = outputs_folder)

# time to threshold
plot_time_to_threshold(outputsfolder = outputs_folder)

# efficacy - saved figure in folder and outputs efficacy results 
hr_results <- plot_hazard_ratios(outputsfolder = outputs_folder)

# plot 1- IRR
plot_irr(outputsfolder = outputs_folder)

# intiail merozoites (only if we export parasitemia)
plot_initial_merozoites(outputsfolder = 'outputs_2025-12-08_treat_0.9start_141threshold5000')
