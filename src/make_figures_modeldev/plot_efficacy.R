# Plot efficacy by year

# Load packages
library(zoo)
library(survival)
library(survminer)
library(broom)
library(ggplot2)

source("R:/Kelly/synergy_orderly/shared/get_cox_efficacy.R")
source("R:/Kelly/synergy_orderly/shared/likelihood.R")

path <- 'R:/Kelly/synergy_orderly/src/sim_cohort_generic/outputs/'
outputsfolder <- 'outputs_2025-11-24'

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

efficacy_plot <- ggplot(tidy_results %>% filter(year == 'Overall')) +
  geom_point(aes(x = term, y = VE, color = reference), size = 2) +
  geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper, color = reference), width = 0.2, linewidth = 1) +
  geom_hline(aes(yintercept = 0)) +
  labs(y = "Vaccine efficacy",
       x = NULL, 
       color = 'Comparison group') +
  scale_color_manual(values = c(
    'SMC reference' = '#2292A4',
    'RTS,S reference' = '#BDBF09',
    'No intervention\nreference' = '#D96C06')) + 
  scale_y_continuous(breaks = seq(min(floor(tidy_results$VE_lower * 10)/10), 1.1, 0.2),
                     limits = c(min(floor(tidy_results$VE_lower * 10)/10),1)) +
  # facet_wrap(~ year) + 
  theme_bw(base_size = 12) 

efficacy_plot

