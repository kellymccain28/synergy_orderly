# Script to run all figure-producing scripts that analyse the outputs of the cohort simulations
setwd('R:/Kelly/synergy_orderly/src/make_figures_modeldev')

filestosource <- list.files(full.names = TRUE)
sapply(filestosource[6:10], source)

outputs_folder <- "outputs_2025-12-02_2_treat_1start_137threshold5000"
outputs_folder <- "outputs_2025-12-02_3_treat_0.9start_150threshold5000"
outputs_folder <- "outputs_2025-12-02_treat_0start_137threshold5000"
outputs_folder <- "outputs_2025-12-02_treat_0.9start_137threshold10000"
outputs_folder <- "outputs_2025-12-05_treat_0.9start_40threshold5000"
# outputs_folder <- "outputs_2025-12-08_treat_0.9start_141threshold5000"
outputs_folder <- "outputs_2025-12-08_treat_0.9start_140threshold5000"

# first, do monthly inci, then order doesn't matter
plot_monthly_incidence(outputsfolder = outputs_folder)

# time to threshold
plot_time_to_threshold(outputsfolder = outputs_folder)

# efficacy - saved figure in folder and outputs efficacy results 
eff_results <- plot_efficacy(outputsfolder = outputs_folder)

# plot 1- IRR
plot_irr(outputsfolder = outputs_folder)

# intiail merozoites (only if we export parasitemia)
plot_initial_merozoites(outputsfolder = 'outputs_2025-12-08_treat_0.9start_141threshold5000')
