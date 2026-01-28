# Reproducing trial outputs from Chandramohan et al 2021
library(survival)
library(survminer)
library(labelled)
library(broom)
library(gtsummary)
library(flexsurv)
library(SurvRegCensCov)
library(haven)
library(cyphr)
library(lubridate)
library(purrr)
library(orderly)
library(MASS)
library(tidyr)
library(tidyverse)

subfolder <- 'figures/'

key <- cyphr::data_key()

orderly_strict_mode()

orderly_dependency(name = 'clean_trial_data',
                   "latest()",
                   files = c('data/primary_persontime.rds',
                             'data/children.rds',
                             'data/mitt.rds',
                             'data/delivery_detail.rds'))

orderly_artefact(files = 'surv_analysis_trial.rds')

orderly_shared_resource('get_cox_efficacy.R')
source('get_cox_efficacy.R')

orderly_shared_resource('get_incidence.R')
source('get_incidence.R')

# Read in the data ----
primary_pt <- cyphr::decrypt(readRDS('data/primary_persontime.rds'), key)
children <- cyphr::decrypt(readRDS('data/children.rds'), key)
mitt <- cyphr::decrypt(readRDS('data/mitt.rds'), key) %>%
  labelled:::remove_labels(mitt)
delivery <- cyphr::decrypt(readRDS('data/delivery_detail.rds'), key)

primary <- primary_pt %>%
  left_join(delivery)

# Survival analysis to reproduce results from trial
## SMC comparator by year and overall ---- 
smcrefresults <- get_cox_efficacy(df = primary_pt, 
                                  ref = 'arm_smcref',
                                  model = FALSE)

## RTSS comparator by year and overall ---- 
rtssrefresults <- get_cox_efficacy(df = primary_pt, 
                                   ref = 'arm_rtssref',
                                   model = FALSE)

# Plot the vaccine efficacies 
tidy_results <- rbind(smcrefresults, rtssrefresults) %>%
  mutate(year = factor(year, 
                levels = c(1, 2, 3, 'overall'),
                labels = c("Year 1", "Year 2", "Year 3", "Overall")))

saveRDS(tidy_results, 'surv_analysis_trial.rds')

# Plot efficacy
efficacies <- ggplot(tidy_results %>% 
                       filter(term %in% c("Both vs. RTSS",'RTSS vs. SMC','Both vs. SMC')))+
  geom_point(aes(x = term, y = VE, group = year, color = year),
             position = position_dodge(width=0.3), size = 2) + 
  geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper, group = year, color = year), 
                position = position_dodge(width=0.3), width = 0.2) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_y_continuous(breaks = seq(-0.4, 1, 0.2),
                     limits = c(-0.3, 1),
                     labels = scales::percent) +
  scale_color_manual(values = c('Year 1' = '#7FB800',
                                'Year 2' = '#573280',
                                'Year 3' = '#CE6479',
                                'Overall' = '#197BBD')) +
  labs(x = '',
       y = 'Efficacy',
       color = 'Trial year') +
  theme_bw(base_size = 14) + 
  theme(legend.position = c(0.85, 0.8),
        legend.background = element_rect(fill = "transparent", color = NA))
efficacies
ggsave(filename = 'efficacy_trial.png', efficacies, height = 8, width = 8)


# Get efficacy for 3 versus 2 doses 
df <- primary %>%
  dplyr::select(-tidyselect::starts_with('y')) %>%
  rowwise() %>%
  mutate(year = lubridate::year(dcontact),
         n_vaccine_doses = sum(as.numeric(nprimary), as.numeric(boost1_done), as.numeric(boost2_done), na.rm = TRUE),
         n_vaccine_doses = ifelse(n_vaccine_doses <= 2, '2 or less', as.character(n_vaccine_doses))) %>%
  mutate(nprimary = factor(nprimary, levels = c(3, 2, 1))) %>%
  mutate(n_vaccine_doses = factor(n_vaccine_doses, levels = c(5, 4, 3, '2 or less'))) %>%
  filter(arm == 'both'| arm == 'rtss') 

coxdoses_adjustsmc <- coxph(Surv(start_time, end_time, event) ~ factor(n_vaccine_doses) + factor(arm) + country, #+ strata(country),# + factor(country),
                  data = df ,
                  cluster = rid,
                  ties = "efron")
coxdoses_interactsmc <- coxph(Surv(start_time, end_time, event) ~ factor(nprimary) * arm + strata(country),
                            data = df,
                            cluster = rid,
                            ties = "efron")

results <- tidy(coxdoses_interactsmc,
                exponentiate = TRUE,
                conf.int = TRUE) %>%
  # Calculate vaccine efficacy and its confidence intervals
  mutate(VE = (1 - estimate) * 100,              # VE = 1 - HR
         VE_lower = (1 - conf.high) * 100,       # Lower CI for VE = 1 - Upper CI for HR
         VE_upper = (1 - conf.low) * 100) %>%     # Upper CI for VE = 1 - Lower CI for HR
  mutate(n_events = coxdoses_interactsmc$nevent,
         n_obs = coxdoses_interactsmc$n)

ggplot(results )+#%>% filter(term != 'factor(arm)rtss'))+
  geom_point(aes(x = term, y = VE)) + 
  geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper), width = 0.2) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  labs(x = '') +
  theme_bw(base_size = 16) 

# Stratified analysis
results_stratified <- df %>%
  # filter(country == 'Mali') %>%
  # filter(dcontact < as.Date('2018-04-01')) %>%
  group_by(arm) %>%
  do({
    mod <- coxph(
      Surv(start_time, end_time, event) ~ n_vaccine_doses,#nprimary,
      data = .,
      cluster = rid,
      ties = "efron"
    )
    tidy(mod, exponentiate = TRUE, conf.int = TRUE)
  }) %>%
  # filter(grepl("nprimary", term)) %>%
  mutate(
    VE = (1 - estimate),
    VE_lower = (1 - conf.high),
    VE_upper = (1 - conf.low),
    doses = case_when(
      term == "nprimary2" ~ "2 doses",
      term == "nprimary1" ~ "1 dose",
      term == 'n_vaccine_doses4' ~ "4 doses",
      term == 'n_vaccine_doses3' ~ "3 doses",
      term == 'n_vaccine_doses2 or less' ~ '2 or fewer doses'
    ),
    arm_label = ifelse(arm == "rtss", "RTS,S only", "RTS,S + SMC")
  )

# Stratified plot
ggplot(results_stratified, aes(x = doses, y = VE, color = arm_label)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) + 
  geom_errorbar(aes(ymin = VE_lower, ymax = VE_upper), 
                width = 0.2, position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
  # scale_x_discrete(limits = c("2 doses", "1 dose")) +
  scale_color_manual(values = c("RTS,S only" = "#E69F00", 
                                "RTS,S + SMC" = "#56B4E9")) +
  scale_y_continuous(labels = scales::percent_format(),
                     breaks = seq(-4, 1, 0.5)) +
  labs(
    x = 'Number of RTS,S doses',
    y = 'VE (%)',
    color = 'Study Arm'
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0),
    legend.position = "bottom"
  )


# The estimate for 2 doses compared to 1 dose is -80%, so 1 dose seems to be more (?) efficacious 
# For 3 doses, it is even worse, and the confidence interval doesn't overlap 0 . --- these first two don't matter 
# Both combination and rtss groups cross 0 which indicates not a large difference with SMC, but rtss seems worse 
# For RTSS+SMC versus SMC, 2 doses is 62.1% better and 3 doses is 63.8% better than 1 dose. 
# For RTSS versus SMC, 2 doses is 11.5% better and 3 doses is 23.5% better than 1 dose, but both cross 0
# dose doesn't seem to matter very much.


# Survival curve
kmsurvobj <- survfit(Surv(start_time, 
                          end_time, 
                          event) ~ arm ,#+ strata(country), 
                     data = primary_pt)

km_plot <- ggsurvplot(kmsurvobj, #group.by = "country",
           # linetype = "strata",
           facet.by = 'country',
           tables.theme = theme_cleantable(),
           conf.int = TRUE,
           fun = 'event',
           risk.table = TRUE,
           cumevents = TRUE,
           ggtheme = theme_bw(),
           palette = 'Dark2')


ggsave(filename = 'km_trial.png', km_plot)


# Get monthly incidence 
monthly_inci <- get_incidence(model = FALSE, 
                              df_children = children, 
                              casedata = mitt)
saveRDS(monthly_inci, 'monthly_incidence_trial.rds')


monthlyincidenceplot <- monthly_inci %>%
  ggplot(aes(x = date, y = incidence_per_1000pm)) +
  geom_point(aes(color = arm), size = 2) +
  geom_ribbon(aes(ymin = lower_per_1000, ymax = upper_per_1000, fill = arm),
                alpha = 0.5) + #, width = 5, linewidth = 1
  geom_line(aes(color = arm), linewidth = 1) +
  # facet_wrap(~arm, nrow = 3) +
  scale_color_manual(values =  c('both' = '#E15554', 
                                 'none' = '#E1BC29',
                                 'rtss' = '#3BB273',
                                 'smc' = '#7768AE',
                                 'SMC delivery' = '#4D9DE0',
                                 'RTS,S delivery' = '#470024'))+#c('#C44536','#772E25','#197278','#283D3B'))+
  scale_fill_manual(values =  c('both' = '#E15554', 
                                 'none' = '#E1BC29',
                                 'rtss' = '#3BB273',
                                 'smc' = '#7768AE'))+#c('#C44536','#772E25','#197278','#283D3B'))+
  scale_y_continuous(breaks = seq(0,150,25)) +
  scale_x_date(breaks = '2 months',
               labels = scales::label_date_short()) + 
  labs(
    # title = "Monthly malaria incidence per 1000 person-months",
    x = "Month",
    y = "Incidence per 1000 person-months",
    color = "Intervention Arm",
    fill = "Intervention Arm"
  ) +
  theme_minimal(base_size = 14)
monthlyincidenceplot
ggsave("trial_monthlyincidence.png", plot = monthlyincidenceplot, bg = 'white', width = 12, height = 6)
ggsave("trial_monthlyincidence.pdf", plot = monthlyincidenceplot, width = 12, height = 6)

# Average delivery times for each intervention / country 
ggplot(delivery %>% filter(arm !='rtss')) +
  geom_histogram(aes(x = nsmc_received, group = arm, fill = arm),
                 position = 'dodge') + 
  facet_wrap(~ country) +
  scale_x_continuous(breaks =seq(0,12)) + 
  labs(x = 'Number of SMC doses received (maximum of 12)',
       y = 'Number of people',
       fill = 'Intervention arm') + 
  theme_bw()


# percent of people in each arm who missed at least 1 round of smc 
delivery %>% 
  filter(arm != 'rtss') %>%
  mutate(missedanysmc = ifelse(nsmc_received == 12, FALSE, TRUE),
         missed3ormoresmc = ifelse(nsmc_received > 9, FALSE, TRUE)) %>%
  group_by(arm, country) %>%
  summarize(missedanysmc = sum(missedanysmc),
            missed3ormoresmc = sum(missed3ormoresmc),
            total = n()) %>% 
  mutate(percent_missinganysmc = scales::percent(missedanysmc / total),
         percent_missed3ormoresmc = scales::percent(missed3ormoresmc / total)) %>%
  arrange(country)
# in Mali, about 41-42% of people missed at least 1 smc round over 3 years in both SMC and Both groups,
# Mali: about 19-20% missed at least 3 rounds
# while in BF, about 41% in the Both group missed at least 1 smc round, while 61% in teh SMC only group missed at least 1 round 
# BF: 15.6% missed at least 3 rounds in both group and 23.1% missed rounds in SMC only group 
y1cols <- grep("^y1", names(delivery), value = TRUE)
y1cols <- y1cols[!(y1cols %in% grep('date', y1cols, value = TRUE))]
delivery$nsmc_doses <- rowSums(!is.na(delivery[, y1cols]))

# Regression to see if the number of SMC doses is related to the arm
mm <- glm(formula = nsmc_doses ~ arm,
          family = gaussian, 
          data = delivery %>% filter(arm != 'rtss', country == 'BF'))
tidy(mm, exponentiate = TRUE, conf.int = TRUE)

# Calculate variance-to-mean ratio by group
delivery %>%
  group_by(country, arm) %>%
  summarise(
    mean = mean(nsmc_doses, na.rm = TRUE),
    variance = var(nsmc_doses, na.rm = TRUE),
    ratio = variance / mean,
    n = n()
  )
# because BF is overdispersed (ratio>>1), I should use the NB distribution instead of Poisson 

glm_stratified <- delivery %>%
  filter(arm != 'rtss') %>%
  group_by(country) %>%
  do({
    # mod <- glm(
    #   nsmc_doses ~ arm,
    #   family = poisson(link = "log"),
    #   data = .,
    #   control = glm.control(maxit = 100)
    # )
    # tidy(mod, exponentiate = TRUE, conf.int = TRUE)
    mod_nb <- glm.nb(nsmc_doses ~ arm, 
                     data = .,
                     control = glm.control(maxit = 50))
    tidy(mod_nb, exponentiate = TRUE, conf.int = TRUE)
  }) %>%
  filter(term == "armsmc") %>%
  mutate(pct_reduction = (1 - estimate) * 100) %>%
  dplyr::select(country, estimate, conf.low, conf.high, pct_reduction, p.value)
glm_stratified

# Visualization
ggplot(glm_stratified, 
       aes(x = country, y = estimate, color = country)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_y_continuous(limits = c(0.9, 1.05), breaks = seq(0.9, 1.05, 0.05)) +
  labs(
    x = NULL,
    y = "Rate Ratio (SMC v both)",
    # title = "SMC Dose Receipt: SMC-only vs Combined Intervention",
    # subtitle = "RR < 1 indicates fewer doses in SMC-only group",
    caption = "Reference group: both\n RR < 1 indicates fewer doses in SMC group."
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

delivery_avg <- delivery %>%
  group_by(country) %>%
  summarize(across(contains('date'),
                   list(median = ~median(.x, na.rm = TRUE)),
                   .names = "{.col}_{.fn}")) #%>%
  # pivot_longer(cols = v1_date_median:fu_end_date_median)
vax_dates_avg <- delivery_avg %>%
  dplyr::select(country, v1_date_median, v2_date_median, v3_date_median, boost1_date_median, boost2_date_median) %>%
  pivot_longer(cols = c(v1_date_median, v2_date_median, v3_date_median, boost1_date_median, boost2_date_median),
               names_to = 'dose',
               values_to = 'date') %>%
  mutate(dose = factor(str_replace(dose, "_date_median", ""), levels = c("v1",'v2','v3','boost1', 'boost2')))
saveRDS(vax_dates_avg, 'R:/Kelly/synergy_orderly/shared/median_rtss_dates.rds')


vaxdates <- delivery %>%
  pivot_longer(cols = c(v1_date, v2_date, v3_date, boost1_date, boost2_date),
               names_to = 'dose',
               values_to = 'date') %>%
  mutate(dose = factor(str_replace(dose, "_date", ""), levels = c("v1",'v2','v3','boost1', 'boost2')),
         date = ifelse(country == 'BF' & dose == 'boost2', 1400, date))

ggplot(delivery_avg) +
  geom_histogram(data = vaxdates, aes(x = date, fill = dose)) +
  scale_x_date(date_breaks = '1 months', labels = scales::label_date_short()) +
  geom_vline(data = delivery_avg, aes(xintercept = v1_date_median, color = 'v1'))+
  geom_vline(data = delivery_avg, aes(xintercept = v2_date_median, color = 'v2'))+
  geom_vline(data = delivery_avg, aes(xintercept = v3_date_median, color = 'v3')) +
  geom_vline(data = delivery_avg, aes(xintercept = boost1_date_median, color = 'boost1')) +
  geom_vline(data = delivery_avg, aes(xintercept = boost2_date_median, color = 'boost2')) +
  # theme(axis.text.x = element_text(angle = 45)) +
  facet_wrap(~country)

median(delivery$y1p1d1_date_received, na.rm = TRUE) # 7-27
median(delivery$y1p2d1_date_received, na.rm = TRUE) # 8-24 
median(delivery$y1p3d1_date_received, na.rm = TRUE) # 9-23

smcdates <- delivery_avg %>%
  select(country, contains('d3')) %>%
  pivot_longer(cols = contains('date_received'),
               names_to = 'smcdose',
               values_to = 'date') 
ggplot(smcdates) +
  geom_vline(aes(xintercept = date, color = smcdose)) + facet_wrap(~country)

saveRDS(smcdates, 'R:/Kelly/synergy_orderly/shared/median_smc_dates.rds')








