# Script to run all figure-producing scripts that analyse the outputs of the cohort simulations
setwd('R:/Kelly/synergy_orderly/src/make_figures_modeldev')

outputs_folder <- 'outputs_2025-12-01_2'

filestosource <- list.files(full.names = TRUE)
sapply(filestosource[6:10], source)

# first, do monthly inci, then order doesn't matter
plot_monthly_incidence(outputsfolder = outputs_folder)

# intiail merozoites (only if we export parasitemia)
plot_initial_merozoites(outputsfolder = 'outputs_2025-11-26')

# time to threshold
plot_time_to_threshold(outputsfolder = outputs_folder)

# efficacy
plot_efficacy(outputsfolder = outputs_folder)

# plot 1- IRR
plot_irr(outputsfolder = outputs_folder)
