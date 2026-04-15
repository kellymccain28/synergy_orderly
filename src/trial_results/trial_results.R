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
library(janitor)
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
                             'data/delivery_detail.rds',
                             'data/weekly.rds',
                             'data/serology.rds'))

orderly_artefact(files = 'surv_analysis_trial.rds')

orderly_shared_resource('get_cox_efficacy.R')
source('get_cox_efficacy.R')

orderly_shared_resource('get_incidence.R')
source('get_incidence.R')

mycols <- c('both' = '#E15554', 
            'none' = '#E1BC29',
            'rtss' = '#3BB273',
            'smc' = '#7768AE',
            'SMC delivery' = '#4D9DE0',
            'RTS,S delivery' = '#470024')

# Read in the data ----
primary_pt <- cyphr::decrypt(readRDS('data/primary_persontime.rds'), key)
children <- cyphr::decrypt(readRDS('data/children.rds'), key)
mitt <- cyphr::decrypt(readRDS('data/mitt.rds'), key) %>%
  labelled:::remove_labels(mitt)
delivery <- cyphr::decrypt(readRDS('data/delivery_detail.rds'), key)
weekly <- cyphr::decrypt(readRDS('data/weekly.rds'), key)
sero <- cyphr::decrypt(readRDS('data/serology.rds'), key)

primary <- primary_pt %>%
  left_join(delivery)

# Survival analysis to reproduce results from trial ----
## SMC comparator by year and overall ---- 
smcrefresults <- get_cox_efficacy(df = primary_pt, 
                                  ref = 'arm_smcref',
                                  model = FALSE)

## RTSS comparator by year and overall ---- 
rtssrefresults <- get_cox_efficacy(df = primary_pt, 
                                   ref = 'arm_rtssref',
                                   model = FALSE)

# Plot the vaccine efficacies ----
tidy_results <- rbind(smcrefresults, rtssrefresults) %>%
  mutate(year = factor(year, 
                       levels = c(1, 2, 3, 'overall'),
                       labels = c("Year 1", "Year 2", "Year 3", "Overall")))

saveRDS(tidy_results, 'surv_analysis_trial.rds')

# Survival analysis to reproduce results from trial by country----
## SMC comparator by year and overall ---- 
BFsmcrefresults <- get_cox_efficacy(df = primary_pt %>% filter(country == 'BF'), 
                                  ref = 'arm_smcref',
                                  model = FALSE)

## RTSS comparator by year and overall ---- 
BFrtssrefresults <- get_cox_efficacy(df = primary_pt %>% filter(country == 'BF'), 
                                   ref = 'arm_rtssref',
                                   model = FALSE)

## SMC comparator by year and overall ---- 
Malismcrefresults <- get_cox_efficacy(df = primary_pt %>% filter(country == 'Mali'), 
                                    ref = 'arm_smcref',
                                    model = FALSE)

## RTSS comparator by year and overall ---- 
Malirtssrefresults <- get_cox_efficacy(df = primary_pt %>% filter(country == 'Mali'), 
                                     ref = 'arm_rtssref',
                                     model = FALSE)

mali_tidy <- rbind(Malirtssrefresults, Malismcrefresults) %>%
  mutate(country= 'Mali')
bf_tidy <- rbind(BFrtssrefresults, BFsmcrefresults)%>%
  mutate(country= 'BF')

tidy_results_stratified <- rbind(mali_tidy, bf_tidy) %>%
  mutate(year = factor(year, 
                       levels = c(1, 2, 3, 'overall'),
                       labels = c("Year 1", "Year 2", "Year 3", "Overall")))

saveRDS(tidy_results_stratified, 'surv_analysis_trial_stratified.rds')

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
  theme(legend.position = c(0.85, 0.75),
        legend.background = element_rect(fill = "transparent", color = NA))
# efficacies
ggsave(filename = 'efficacy_trial.pdf', efficacies, height = 5, width = 5)


efficacies_strat <- ggplot(tidy_results_stratified %>% 
                       filter(term %in% c("Both vs. RTSS",'RTSS vs. SMC','Both vs. SMC')))+
  geom_point(aes(x = term, y = VE, group = year, color = year),
             position = position_dodge(width=0.3), size = 2) + 
  geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper, group = year, color = year), 
                position = position_dodge(width=0.3), width = 0.2) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_y_continuous(breaks = seq(-0.4, 1, 0.2),
                     limits = c(-0.4, 1),
                     labels = scales::percent) +
  scale_color_manual(values = c('Year 1' = '#7FB800',
                                'Year 2' = '#573280',
                                'Year 3' = '#CE6479',
                                'Overall' = '#197BBD')) +
  facet_wrap(~country) +
  labs(x = '',
       y = 'Efficacy',
       color = 'Trial year') +
  theme_bw(base_size = 14) + 
  theme(legend.position = c(0.9, 0.8),
        legend.background = element_rect(fill = "transparent", color = NA))
# efficacies_strat
ggsave(filename = 'efficacy_trial_stratified.pdf', efficacies_strat, height = 5, width = 10)


# Get efficacy for 3 versus 2 doses 
df <- primary %>%
  dplyr::select(-tidyselect::starts_with('y')) %>%
  rowwise() %>%
  mutate(year = lubridate::year(dcontact),
         n_vaccine_doses = sum(as.numeric(nprimary), as.numeric(boost1_done), as.numeric(boost2_done), na.rm = TRUE),
         n_vaccine_doses = ifelse(n_vaccine_doses <= 2, '2 or fewer', as.character(n_vaccine_doses))) %>%
  mutate(nprimary = factor(nprimary, levels = c(3, 2, 1))) %>%
  mutate(n_vaccine_doses = factor(n_vaccine_doses, levels = c(5, 4, 3, '2 or fewer'))) %>%
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

# ggplot(results )+#%>% filter(term != 'factor(arm)rtss'))+
#   geom_point(aes(x = term, y = VE)) + 
#   geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper), width = 0.2) +
#   geom_hline(aes(yintercept = 0), linetype = 2) +
#   labs(x = '') +
#   theme_bw(base_size = 16) 

# Stratified analysis
results_stratified <- df %>%
  # filter(country == 'BF') %>%
  # filter(dcontact < as.Date('2018-04-01')) %>%
  group_by(arm) %>%
  do({
    mod <- coxph(
      Surv(start_time, end_time, event) ~ nprimary + strata(country),#n_vaccine_doses,#
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
      term == 'n_vaccine_doses2 or fewer' ~ '2 or fewer doses'
    ),
    arm_label = ifelse(arm == "rtss", "RTS,S only", "RTS,S + SMC")
  )

# Stratified plot
ve_strat <- ggplot(results_stratified, aes(x = doses, y = VE, color = arm_label)) +
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


# Get monthly incidence ----
monthly_inci <- get_incidence(model = FALSE, 
                              df_children = children, 
                              casedata = mitt)
saveRDS(monthly_inci, 'monthly_incidence_trial.rds')

monthly_inci_BF <- get_incidence(model = FALSE, 
                                 df_children = children %>% filter(country == 'BF'), 
                                 casedata = mitt %>% filter(country == 'BF'))
saveRDS(monthly_inci_BF, 'monthly_incidence_trial_BF.rds')

monthly_inci_Mali <- get_incidence(model = FALSE, 
                                   df_children = children %>% filter(country == 'Mali'), 
                                   casedata = mitt %>% filter(country == 'Mali'))
saveRDS(monthly_inci_Mali, 'monthly_incidence_trial_Mali.rds')


# Get median smc and rtss dates ----
smc_dates <- readRDS('R:/Kelly/synergy_orderly/shared/median_smc_dates.rds') 
smc_lines_BF <- smc_dates %>% ungroup() %>%
  filter(country == 'BF' & arm != 'rtss') %>%
  dplyr::select(date, arm) %>% 
  mutate(color = '#709176') 
smc_lines_Mali <-  smc_dates %>% ungroup() %>%
  filter(country == 'Mali' & arm != 'rtss') %>%
  dplyr::select(date, arm) %>% 
  mutate(color = '#709176') 

vaxdates <- readRDS('R:/Kelly/synergy_orderly/shared/median_rtss_dates.rds')
rtss_lines_BF <-  vaxdates %>% ungroup() %>%
  filter(country == 'BF' & arm != 'smc') %>%
  dplyr::select(date, arm) %>% 
  mutate(color = '#59114D') 
rtss_lines_Mali <-  vaxdates %>% ungroup() %>%
  filter(country == 'Mali' & arm != 'smc') %>%
  dplyr::select(date, arm) %>% 
  mutate(color = '#59114D') 


# make inci plots ----
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
# monthlyincidenceplot
ggsave("trial_monthlyincidence.png", plot = monthlyincidenceplot, bg = 'white', width = 12, height = 6)
ggsave("trial_monthlyincidence.pdf", plot = monthlyincidenceplot, width = 12, height = 6)

#Burkina only 
monthlyincidenceplotbf <- monthly_inci_BF %>%
  ggplot(aes(x = date, y = incidence_per_1000pm)) +
  geom_point(aes(color = arm), size = 2) +
  geom_ribbon(aes(ymin = lower_per_1000, ymax = upper_per_1000, fill = arm),
              alpha = 0.5) + #, width = 5, linewidth = 1
  geom_line(aes(color = arm), linewidth = 1) +
  geom_vline(data = smc_lines_BF, aes(xintercept = date, color = 'SMC delivery', group = arm), linetype = 2, linewidth = 0.7)+
  geom_vline(data = rtss_lines_BF, aes(xintercept = date, color = 'RTS,S delivery', group = arm), linetype = 3, linewidth = 0.8) +
  # facet_wrap(~arm, nrow = 3) +
  scale_color_manual(values =  c('both' = '#E15554', 
                                 'none' = '#E1BC29',
                                 'rtss' = '#3BB273',
                                 'smc' = '#7768AE',
                                 'SMC delivery' = '#4D9DE0',
                                 'RTS,S delivery' = '#470024'),
                     breaks = c('none', 'rtss', 'smc', 'both','SMC delivery','RTS,S delivery'))+#c('#C44536','#772E25','#197278','#283D3B'))+
  scale_fill_manual(values =  c('both' = '#E15554', 
                                'none' = '#E1BC29',
                                'rtss' = '#3BB273',
                                'smc' = '#7768AE'),
                    guide = 'none')+#c('#C44536','#772E25','#197278','#283D3B'))+
  scale_y_continuous(breaks = seq(0,max(monthly_inci_BF$upper_per_1000),25)) +
  scale_x_date(breaks = '2 months',
               labels = scales::label_date_short()) + 
  labs(
    # title = "Monthly malaria incidence per 1000 person-months",
    x = "Month",
    y = "Incidence per 1000 person-months",
    color = NULL
  ) +
  theme_minimal(base_size = 14) + 
  facet_wrap(~factor(arm, levels = c('none','rtss','smc','both')),
             nrow = 3)

ggsave("trial_monthlyincidence_BF.png", plot = monthlyincidenceplotbf, bg = 'white', width = 12, height = 8)
ggsave("trial_monthlyincidence_BF.pdf", plot = monthlyincidenceplotbf, width = 12, height = 8)

#Mali only 
monthlyincidenceplotmali <- monthly_inci_Mali %>%
  ggplot(aes(x = date, y = incidence_per_1000pm)) +
  geom_point(aes(color = arm), size = 2) +
  geom_ribbon(aes(ymin = lower_per_1000, ymax = upper_per_1000, fill = arm),
              alpha = 0.5) + #, width = 5, linewidth = 1
  geom_line(aes(color = arm), linewidth = 1) +
  geom_vline(data = smc_lines_BF, aes(xintercept = date, color = 'SMC delivery', group = arm), linetype = 2, linewidth = 0.7)+
  geom_vline(data = rtss_lines_BF, aes(xintercept = date, color = 'RTS,S delivery', group = arm), linetype = 3, linewidth = 0.8) +
  # facet_wrap(~arm, nrow = 3) +
  scale_color_manual(values =  c('both' = '#E15554', 
                                 'none' = '#E1BC29',
                                 'rtss' = '#3BB273',
                                 'smc' = '#7768AE',
                                 'SMC delivery' = '#4D9DE0',
                                 'RTS,S delivery' = '#470024'),
                     breaks = c('none', 'rtss', 'smc', 'both','SMC delivery','RTS,S delivery'))+#c('#C44536','#772E25','#197278','#283D3B'))+
  scale_fill_manual(values =  c('both' = '#E15554', 
                                'none' = '#E1BC29',
                                'rtss' = '#3BB273',
                                'smc' = '#7768AE',
                                'SMC delivery' = '#4D9DE0',
                                'RTS,S delivery' = '#470024'),
                    guide = "none" )+#c('#C44536','#772E25','#197278','#283D3B'))+
  scale_y_continuous(breaks = seq(0,max(monthly_inci_Mali$upper_per_1000),25)) +
  scale_x_date(breaks = '2 months',
               labels = scales::label_date_short()) + 
  labs(
    # title = "Monthly malaria incidence per 1000 person-months",
    x = "Month",
    y = "Incidence per 1000 person-months",
    color = NULL#"Intervention Arm"
  ) +
  theme_minimal(base_size = 14) + 
  facet_wrap(~factor(arm, levels = c('none','rtss','smc','both')),
             nrow = 3)
ggsave("trial_monthlyincidence_Mali.png", plot = monthlyincidenceplotmali, bg = 'white', width = 12, height = 8)
ggsave("trial_monthlyincidence_Mali.pdf", plot = monthlyincidenceplotmali, width = 12, height = 8)


# Delivery proportion ----
nsmc_percountry <- delivery %>%
  ungroup() %>%
  group_by(arm, country, nsmc_received) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(country, arm) %>%
  mutate(psmc_received = count / sum(count))

#smc 
nsmcplot <- ggplot(nsmc_percountry %>% filter(arm !='rtss')) +
  geom_col(aes(x = as.factor(nsmc_received), y = psmc_received, 
                     group = arm, fill = arm),
                 position = 'dodge') + 
  facet_wrap(~ country) +
  # scale_x_continuous(breaks =seq(0,12)) + 
  scale_y_continuous(labels = scales::percent,
                     breaks= seq(0,0.6,.10)) +
  scale_fill_manual(values = mycols) +
  labs(x = 'Number of SMC rounds received (maximum of 12)',
       y = 'Percent of children',
       fill = 'Intervention arm') + 
  theme_bw(base_size =  14)
ggsave('nsmc_received_bycountryandarm.pdf', plot = nsmcplot, height = 5, width = 10)

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
         percent_missed3ormoresmc = scales::percent(missed3ormoresmc / total),
         atleast11 = scales::percent(1-(missedanysmc / total)),
         atleast9 = scales::percent(1-(missed3ormoresmc / total))) %>%
  arrange(country)
# in Mali, about 41-42% of people missed at least 1 smc round over 3 years in both SMC and Both groups,
# Mali: about 19-20% missed at least 3 rounds
# while in BF, about 41% in the Both group missed at least 1 smc round, while 61% in teh SMC only group missed at least 1 round 
# BF: 15.6% missed at least 3 rounds in both group and 23.1% missed rounds in SMC only group 
delivery %>%
  filter(arm!= 'rtss') %>%
  tabyl(nsmc_received, country) %>%
  adorn_percentages() %>%
  adorn_pct_formatting()

ycols <- grep("^y", names(delivery), value = TRUE)
ycols <- ycols[!(ycols %in% grep('date', ycols, value = TRUE))]
delivery$nsmc_doses <- rowSums(delivery[, ycols], na.rm = TRUE)

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
rr_doses <- ggplot(glm_stratified, 
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


#rtss  
delivery <- delivery %>%
  rowwise() %>%
  mutate(n_rtss =  sum(as.numeric(nprimary), as.numeric(boost1_done), as.numeric(boost2_done), na.rm = TRUE))
nrtss_percountry <- delivery %>%
  ungroup() %>%
  group_by(arm, country, n_rtss) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(country, arm) %>%
  mutate(prtss_received = count / sum(count))

delivery %>%
  filter(arm!= 'smc') %>%
  tabyl(arm, nprimary, country) %>%
  adorn_percentages() %>%
  adorn_pct_formatting()

delivery %>%
  filter(arm!= 'smc') %>%
  tabyl(arm, n_rtss, country) %>%
  adorn_percentages() %>%
  adorn_pct_formatting()
  
#mean doses per arm 
delivery %>%
  filter(arm != 'smc') %>%
  group_by(country, arm) %>%
  summarise(median_rtss = median(n_rtss),
            median_smcdoses = median(nsmc_doses),
            median_smc= median(nsmc_received))
#rtss 
nrtss <- ggplot(nrtss_percountry %>% filter(arm !='smc')) +
  geom_col(aes(x = as.factor(n_rtss), y = prtss_received, 
               group = arm, fill = arm),
           position = 'dodge') + 
  facet_wrap(~ country) +
  # scale_x_continuous(breaks =seq(0,12)) + 
  scale_y_continuous(labels = scales::percent,
                     breaks= seq(0,0.9,.10),
                     limits = c(0,0.9)) +
  scale_fill_manual(values = mycols) +
  labs(x = 'Number of doses of RTS,S received (maximum of 5)',
       y = 'Percent of children',
       fill = 'Intervention arm') + 
  theme_bw(base_size =  14)
ggsave('nrtss_received_bycountryandarm.pdf', plot = nrtss, height = 5, width = 10)




# Delivery times for each intervention / country ----

delivery_avg <- delivery %>%
  group_by(country, arm) %>%
  summarize(across(contains('date'),
                   list(median = ~median(.x, na.rm = TRUE)),
                   .names = "{.col}_{.fn}")) #%>%
# pivot_longer(cols = v1_date_median:fu_end_date_median)
vax_dates_avg_wide <- delivery_avg %>%
  dplyr::select(country, arm, v1_date_median, v2_date_median, v3_date_median, boost1_date_median, boost2_date_median) 
saveRDS(vax_dates_avg_wide, 'R:/Kelly/synergy_orderly/shared/median_rtss_dates_wide.rds')

vax_dates_avg <- vax_dates_avg_wide %>%
  pivot_longer(cols = c(v1_date_median, v2_date_median, v3_date_median, boost1_date_median, boost2_date_median),
               names_to = 'dose',
               values_to = 'date') %>%
  mutate(dose = factor(str_replace(dose, "_date_median", ""), levels = c("v1",'v2','v3','boost1', 'boost2')))
saveRDS(vax_dates_avg, 'R:/Kelly/synergy_orderly/shared/median_rtss_dates.rds')

vaxdates <- delivery %>%
  pivot_longer(cols = c(v1_date, v2_date, v3_date, boost1_date, boost2_date),
               names_to = 'dose',
               values_to = 'date') %>%
  mutate(dose = factor(str_replace(dose, "_date", ""), levels = c("v1",'v2','v3','boost1', 'boost2')))

smcdates <- delivery %>%
  dplyr::select(rid, arm, country, contains('date_received')) %>%
  pivot_longer(cols = c(y1p1d1_date_received:y3p4d3_date_received),
               names_to = 'smcdose',
               values_to = 'date') %>%
  mutate(# Extract components BEFORE reformatting (easier with original format)
    year = as.numeric(str_extract(smcdose, "(?<=y)\\d+")),
    round = as.numeric(str_extract(smcdose, "(?<=p)\\d+")),
    round = paste0('Year ',year, ', Round ', round),
    smcdose = str_replace(smcdose, "_date_received", ""),
         smcdose = str_replace_all(smcdose, 
                                   "y(\\d+)p(\\d+)d(\\d+)", 
                                   "Year \\1, Round \\2, Dose \\3"))


dosecolors <- c('v1'='#99B3FF',
                'v2'='#7A82AB',
                'v3'='#3C908E',
                'boost1' = '#2DC2BD',
                'boost2' = '#0F5744')
armlabs = c('Both','RTS,S only','SMC only')
names(armlabs) = c('both','rtss','smc')

# rtssdelivery_dates <- ggplot(vaxdates %>% filter(arm != 'smc')) +
#   geom_histogram(aes(x = as.Date(date), fill = dose), binwidth = 1) +
#   scale_x_date(date_breaks = '1 months', labels = scales::label_date_short()) +
#   scale_fill_manual(values = dosecolors,
#                     labels = c('v1' = 'Dose 1',
#                                'v2' = 'Dose 2',
#                                'v3' = 'Dose 3',
#                                'boost1' = 'Booster 1',
#                                'boost2' = 'Booster 2'))+
#   # theme(axis.text.x = element_text(angle = 45)) +
#   facet_grid(vars(country, arm), scales = 'free', labeller = labeller(arm= armlabs)) +
#   labs(x = 'Date',
#        y = 'Number of trial participants',
#        fill = NULL) +
#   theme_bw(base_size = 14)
vaxdates_sml <- vaxdates %>%
  mutate(month = month(date),
         week = week(date), 
         year = year(date)) %>%
  filter(month %in% c(4,5,6,7,8,9)) %>%
  group_by(country, arm, dose, week, year) %>%
  summarise(n = n()) %>%
  mutate(date = ISOweek::ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-1")),
         date_f = factor(format(date, "%b %d %Y")),  # convert to factor
         dose = factor(dose, levels = c('boost2','boost1','v3','v2','v1')))

line_positions <- c(2.5, 1.5, 0.5)
rtssdelivery_dates <- ggplot(vaxdates_sml %>% filter(arm != 'smc')) +
  geom_tile(aes(x = date, y = dose, fill = n), color = 'black') +
  geom_hline(yintercept = line_positions, color = "grey50", size = 0.2, linetype = "dashed") +
  scale_x_date(date_breaks = '1 month', 
               labels = scales::label_date_short(),
               expand = c(0,0)) +
  scale_y_discrete(labels = c('v1'= 'Dose 1',
                              'v2'= 'Dose 2',
                              'v3' = 'Dose 3',
                              'boost1' = 'Booster 1',
                              'boost2'='Booster 2'))+
  scale_fill_viridis_c() +
  facet_grid(rows = vars(country, arm), cols = vars(year), 
             scales = 'free_x', labeller = labeller(arm = armlabs)) +
  labs(x = 'Date', y = NULL, fill = 'Number of participants') +
  theme_classic(base_size = 14)
ggsave('rtss_delivery_dates.pdf', plot = rtssdelivery_dates, height = 7, width = 12)

vaxdates2 <- vaxdates %>%
  dplyr::select(rid, arm, dose, date) %>%
  group_by(date, dose, arm) %>%
  count()
# ggplot(vaxdates2 %>% filter(date < '2017-12-01')) + 
#   geom_tile(aes(x = as.Date(date), y = dose, fill = n)) + 
#   theme_classic()

# smcdelivery_dates <- ggplot(smcdates %>% filter(arm != 'rtss')) +
#   geom_histogram(aes(x = as.Date(date), fill = round), binwidth = 1) +
#   scale_x_date(date_breaks = '1 months', labels = scales::label_date_short()) +
#   facet_grid(vars(country, arm), scales = 'free', labeller = labeller(arm= armlabs)) +
#   labs(x = 'Date',
#        y = 'Number of trial participants',
#        fill = NULL) +
#   theme_bw(base_size = 14)
smcdates_sml <- smcdates %>%
  filter(grepl('Dose 1',smcdose)) %>%
  mutate(month = month(date),
         week = week(date), 
         year = year(date)) %>%
  filter(month %in% c(7,8,9,10,11)) %>%
  group_by(country, arm, round, week, year) %>%
  summarise(n = n()) %>%
  mutate(date = ISOweek::ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-1")),
         round= factor(round, levels = c("Year 3, Round 4", "Year 3, Round 3", "Year 3, Round 2", "Year 3, Round 1",
                                              "Year 2, Round 4", "Year 2, Round 3", "Year 2, Round 2", "Year 2, Round 1",
                                              "Year 1, Round 4", "Year 1, Round 3", "Year 1, Round 2", "Year 1, Round 1")))

# Find positions between rounds 
line_positions <- c(8.5, 4.5, 0.5, 8.5, 4.5, 0.5, 8.5, 4.5, 0.5, 8.5, 4.5, 0.5)

smcdelivery_dates <- ggplot(smcdates_sml %>% filter(arm != 'rtss')) +
  geom_tile(aes(x = date, y = round, fill = n), color = 'black') +
  geom_hline(yintercept = line_positions, color = "grey50", size = 0.2, linetype = "dashed") +
  scale_x_date(date_breaks = '1 month', 
               labels = scales::label_date_short(),
               expand = c(0,0)) +
  scale_fill_viridis_c() +
  facet_grid(rows = vars(country, arm), cols = vars(year), 
             scales = 'free_x', labeller = labeller(arm = armlabs)) +
  labs(x = 'Date', y = NULL, fill = 'Number of children') +
  theme_classic(base_size = 14)
ggsave('smcdelivery_dates.pdf', plot = smcdelivery_dates, height = 9, width = 12)


delivery %>% group_by(country) %>% summarise(median(y1p1d1_date_received, na.rm = TRUE)) # 7-27
delivery %>% group_by(country) %>% summarise(median(delivery$y1p2d1_date_received, na.rm = TRUE)) # 8-24 
delivery %>% group_by(country) %>% summarise(median(delivery$y1p3d1_date_received, na.rm = TRUE)) # 9-23
delivery %>% group_by(country) %>% summarise(median(delivery$v1_date, na.rm = TRUE)) 
delivery %>% group_by(country) %>% summarise(median(delivery$v2_date, na.rm = TRUE)) 
delivery %>% group_by(country) %>% summarise(median(delivery$v3_date, na.rm = TRUE)) 

smcdates_wide <- delivery_avg %>%
  dplyr::select(country, arm, contains('d1')) 
saveRDS(smcdates_wide, 'R:/Kelly/synergy_orderly/shared/median_smc_dates_wide.rds')

smcdates <- smcdates_wide %>%
  pivot_longer(cols = contains('date_received'),
               names_to = 'smcdose',
               values_to = 'date') 
# ggplot(smcdates) +
#   geom_vline(aes(xintercept = date, color = smcdose, linetype = arm)) + facet_wrap(~country)

saveRDS(smcdates, 'R:/Kelly/synergy_orderly/shared/median_smc_dates.rds')


# Parasitaemia  ----
weekly <- weekly %>%
  mutate(poutcome= ifelse(pf_asex_fdensity >= 5000, '1', '0'))

parasit_country <- ggplot(weekly %>% filter(!is.na(poutcome))) +
  geom_jitter(aes(x = arm, y = pf_asex_fdensity, color = poutcome)) + 
  geom_violin(aes(x = arm, y= pf_asex_fdensity), alpha = 0.0) +
  geom_hline(yintercept = 5000, linetype = 2) +
  facet_wrap(~country) +
  scale_y_log10() +
  scale_color_manual(values = c('1' = '#F08700',
                                '0' = '#00A6A6'),
                     labels = c('1' = 'Above threshold',
                                '0' = 'Below threshold')) +
  labs(x = NULL,
       y = 'Parasitaemia (PRBCs per \u03bcL)',
       color = NULL) +
  theme_minimal(base_size = 14)
ggsave(filename = 'parasitaemia_bycountry.pdf', plot = parasit_country, height = 4, width = 7)

parasit_overtime <- ggplot(weekly %>% filter(!is.na(pf_asex_fdensity))) +
  geom_point(aes(x = dateweekly, y = pf_asex_fdensity, color = country)) + 
  # geom_smooth(aes(x = dateweekly, y = pf_asex_fdensity, color = country, fill = country), alpha = 0.1) +
  geom_hline(yintercept = 5000, linetype = 2) +
  scale_y_log10() +
  scale_x_date(breaks = '3 months',
               labels = scales::label_date_short()) +
  scale_color_manual(values = c('BF' = '#136F63',
                                'Mali' = '#E0CA3C')) + 
  scale_fill_manual(values = c('BF' = '#136F63',
                                'Mali' = '#E0CA3C')) + 
  facet_wrap(~arm, nrow = 3) + 
  labs(x = 'Date of weekly survey',
       y = 'Parasitaemia (PRBCs per \u03bcL)',
       color = 'Country', fill = 'Country')+
  theme_minimal(base_size = 14)
# year 2, a lot more positive samples it seems
ggsave(filename = 'parasitaemia_overtime.pdf', plot = parasit_overtime, height = 6, width = 8)

weekly %>% 
  group_by(arm) %>%
  summarise(meanpb= mean(pf_asex_fdensity, na.rm = TRUE),
            medianpb = median(pf_asex_fdensity, na.rm = TRUE),
            minpb = min(pf_asex_fdensity, na.rm = TRUE))

# repeated ids
weekly %>% count(rid) 
repeats <- weekly %>%
  filter(rid %in% (weekly %>% count(rid) %>% filter(n > 1) %>% pull(rid))) %>%
  filter(poutcome == 1)
table(repeats$country)
# no ids with repeated parasitaemia results 

# prevalence of parasiteamia in weekly surveys 
weekly %>%
  tabyl(country, poutcome) %>%
  adorn_percentages() %>%
  adorn_pct_formatting()
weekly %>%
  tabyl(poutcome, country) %>%
  adorn_totals() 
# number of non zero samples 
weekly %>%
  tabyl(country,pf_asex_fresult) %>%
  adorn_totals()  %>%
  adorn_percentages() %>%
  adorn_pct_formatting()

weekly %>% 
  group_by(country) %>% 
  summarise(maxdate = max(dateweekly), mindate = min(dateweekly))

# Serology ---- 
ncasesperperson <- mitt %>%
  group_by(rid, country, arm, dcontact) %>%
  count()

sero_withcases <- sero %>%
  left_join(ncasesperperson)

nperyear <- sero %>%
  mutate(year = lubridate::year(sdate)) %>% 
  group_by(arm, country, year) %>% 
  count()

nperpersonperyear <- sero %>%
  mutate(year = lubridate::year(sdate)) %>% 
  group_by(arm, country, rid) %>% 
  count()

summarisedsero <- sero %>% filter(postonly =='post' & arm != 'smc') %>%
  group_by(arm, country, timing) %>%
  summarise(mean = mean(lnnew, na.rm = TRUE),
            median = median(lnnew, na.rm = TRUE))

# ggplot(sero %>% filter(postonly =='post')) + 
#   geom_violin(aes(x = country, y = lnnew, fill = arm))
# 
# ggplot(sero %>% filter(postonly =='post' & arm != 'smc' )) + 
#   geom_point(aes(x = sdate, y = lnnew, color = arm)) + 
#   # geom_line(aes(x = dcontact, y = lnnew, group = rid)) +
#   geom_vline(data = vaxdates, aes(xintercept = date)) + 
#   facet_wrap(~country)


ks_results <- sero %>%
  filter(arm %in% c('both', 'rtss'), postonly == 'post') %>%
  group_by(country, timing) %>%
  summarise(
    # Get values for each arm
    both_values = list(lnnew[arm == 'both']),
    rtss_values = list(lnnew[arm == 'rtss']),
    
    # Sample sizes
    n_both = length(both_values[[1]]),
    n_rtss = length(rtss_values[[1]]),
    
    # Perform KS test
    ks_test = list(ks.test(both_values[[1]], rtss_values[[1]])),
    
    # Extract statistics
    ks_statistic = ks_test[[1]]$statistic,
    ks_pvalue = ks_test[[1]]$p.value,
    
    .groups = 'drop'
  ) %>%
  dplyr::select(-both_values, -rtss_values, -ks_test)  # Remove list columns if desired

# Add p-value annotations to your density plot
ks_results_annot <- ks_results %>%
  mutate(
    label = sprintf("p = %.3f", ks_pvalue),#, ks_statistic),
    max_density = 0.15  # Adjust based on your density plot y-axis
  )

sero_armcountry <- ggplot(sero %>% filter(arm !='smc' & postonly == 'post')) + 
  geom_density(aes(x = lnnew, fill = arm), alpha =0.7) + 
  geom_vline(data = summarisedsero,
             aes(xintercept = median, color = arm), linewidth = 1, linetype = 2, alpha = 0.5) +
  geom_text(data = ks_results_annot,
            aes(x = -Inf, y = Inf, label = label),  # Top left corner
            inherit.aes = FALSE,
            hjust = -0.1,   # Slightly inside from left edge
            vjust = 3,    # Slightly below top edge
            size = 3.5,
            fontface = "italic") +
  scale_fill_manual(values = mycols) + 
  scale_color_manual(values = colorspace::darken(mycols, 0.2)) + 
  facet_wrap(~country + timing) + 
  labs(y = 'Density',
       x = expression(paste(italic('ln'), '(anti-CSP antibody titre in EU/mL)')),
       color = NULL, fill = NULL) +
  theme_minimal(base_size = 14)
ggsave(filename = 'sero_byarmandcountry.pdf', plot = sero_armcountry, height = 5, width = 8)

# Check normality (violated)
shapiro_results <- sero %>%
  filter(arm %in% c('both', 'rtss'), postonly == 'post') %>%
  group_by(country, timing, arm) %>%
  summarise(
    shapiro_p = shapiro.test(lnnew)$p.value,
    normal = shapiro_p > 0.05,
    .groups = 'drop'
  )
print(shapiro_results)  # Most p < 0.05, indicating non-normality

# Create a summary table for reporting
# summary_table <- ks_results %>%
#   mutate(
#     Result = case_when(
#       ks_pvalue > 0.05 ~ "Not significant",
#       ks_pvalue > 0.01 ~ "Significant (p < 0.05)",
#       ks_pvalue > 0.001 ~ "Significant (p < 0.01)",
#       TRUE ~ "Significant (p < 0.001)"
#     ),
#     `KS statistic` = round(ks_statistic, 3),
#     `P-value` = round(ks_pvalue, 4)
#   ) %>%
#   dplyr::select(Country = country, Timing = timing, 
#          `N (both)` = n_both, `N (rtss)` = n_rtss,
#          `KS statistic`, `P-value`, Result)

# Export to CSV
write.csv(summary_table, "ks_test_results_both_vs_rtss.csv", row.names = FALSE)

# linear regression of the measured titre 

sero_wide_v3 <- sero %>% filter(str_detect(timing, 'dose 3')) %>%
  dplyr::select(arm, country, rid, timing, lnnew, age_at_vac) %>%
  pivot_wider(
    names_from = timing,    
    values_from = lnnew     
  ) 
model <- lm(`post dose 3` ~ arm + strata(country) + age_at_vac, data = sero_wide_v3 %>% filter(arm != 'smc'))
summary(model)





# Calculate IRRs as in model outputs ----
source("R:/Kelly/synergy_orderly/src/make_figures_modeldev/bootstrap_metric.R")
library(purrr)


agg_unit = 'halfyear'

calc_irr_trial <- function(agg_unit,
                           armfit,
                           country){
  # should not really matter if i use the bestruns or syntest because i only look at the incidence 
  
  
  # Folders below are from the bestreps runs for the 2-arm fitting (before finished)
  # outputfolder <- c('Mali' = 'outputs_2026-03-26_2', #-- antagonistic
  #                   'BF' = 'outputs_2026-03-26' ## -- synergistic 
  # )
  # Folders below are from the bestreps runs for the 2-arm fitting (after finished)
  if(armfit == 2){ 
    outputfolder <- c('Mali' = 'outputs_2026-03-30_2', #-- antagonistic (slightly)
                      'BF' = 'outputs_2026-03-30' ## -- synergistic
    )}
  # Folders below are from the bestreps runs for the 3-arm fitting 
  if(armfit == 3){ 
    outputfolder <- c('Mali' = 'outputs_2026-03-24_2', # --synergistic
                      'BF' = 'outputs_2026-03-24' # --synergistic - a lot
    )}
  # # Folders below are from the synergy testing runs for the 2-arm fitting(before finished)
  # outputfolder <- c('Mali' = 'outputs_2026-03-26_3', # -- antagonistic
  #                   'BF' = 'outputs_2026-03-26_4' # -- synergistic
  # )
  # # Folders below are from the synergy testing runs for the 2-arm fitting(after finished)
  # outputfolder <- c('Mali' = 'outputs_2026-03-30_8', # -- antagonistic (less than bestreps)
  #                   'BF' = 'outputs_2026-03-30_7' # -- synergistic
  # )
  # # Folders below are from the synergy testing runs for the 3-arm fitting
  # outputfolder <- c('Mali' = 'outputs_2026-03-25_5', # -- synergistic
  #                   'BF' = 'outputs_2026-03-25_6' # --synergistic 
  # )
  outputs_folders <-  c(#'BF3syn' = 'outputs_2026-03-25_6',
    # 'Mali3syn' = 'outputs_2026-03-25_5',
    # 'BF2syn' = 'outputs_2026-03-30_7',
    # 'Mali2syn' = 'outputs_2026-03-30_8',
    
    'BF3best' = 'outputs_2026-03-24',
    'Mali3best' = 'outputs_2026-03-24_2',
    'BF2best' = 'outputs_2026-03-30',
    'Mali2best' = 'outputs_2026-03-30_2')
  
  if(country == 'Mali'){
    # model_inci_summary_all_halfyear <- readRDS(paste0("R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs/",outputfolder['Mali'],"/inci_summary_all_",agg_unit, ".rds")) %>%
    #   mutate(shape_var = case_when(
    #     metric == 'efficacy' & grepl('Expected', comparison) ~ 'Expected',
    #     metric == 'efficacy' ~ 'Model-predicted',
        # TRUE ~ NA)) 
    model_inci_summary_wide_halfyear <- readRDS(paste0("R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs/",outputfolder['Mali'],"/inci_summary_wide_",agg_unit, ".rds"))
    # formatted_infrecords <- readRDS(paste0("R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs/",outputfolder['Mali'],"/formatted_infrecords.rds"))
    
    monthly_incidence_trial <- readRDS("R:/Kelly/synergy_orderly/archive/trial_results/20260324-150302-df3f8aff/monthly_incidence_trial_Mali.rds")
  } else if (country == 'BF'){
    # model_inci_summary_all_halfyear <- readRDS(paste0("R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs/",outputfolder['BF'],"/inci_summary_all_",agg_unit, ".rds")) %>%
    #   mutate(shape_var = case_when(
    #     metric == 'efficacy' & grepl('Expected', comparison) ~ 'Expected',
    #     metric == 'efficacy' ~ 'Model-predicted',
        # TRUE ~ NA)) 
    model_inci_summary_wide_halfyear <- readRDS(paste0("R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs/",outputfolder['BF'],"/inci_summary_wide_",agg_unit, ".rds"))
    # formatted_infrecords <- readRDS(paste0("R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs/",outputfolder['BF'],"/formatted_infrecords.rds"))
    
    monthly_incidence_trial <- readRDS("R:/Kelly/synergy_orderly/archive/trial_results/20260324-150302-df3f8aff/monthly_incidence_trial_BF.rds")
  }
  
  monthly_incidence_trial <- monthly_incidence_trial %>%
    filter(!is.na(date)) %>%
    filter(date > '2017-05-01') 
  
  if(agg_unit == 'halfyear'){
    inci_annual <- monthly_incidence_trial %>%
      mutate(time_value = case_when(date < '2017-10-01' & date > '2017-05-01' ~ 'June 2017-Sep 2017',
                                    date < '2018-04-01' ~ 'Oct 2017-March 2018',
                                    date < '2018-10-01' ~ 'April 2018-Sep 2018',
                                    date < '2019-04-01' ~ 'Oct 2018-March 2019',
                                    date < '2019-10-01' ~ 'April 2019-Sep 2019',
                                    date < '2020-04-01' ~ 'Oct 2019-March 2020'),
             time_value_num = case_when(date < '2017-10-01' & date > '2017-05-01' ~ '1',
                                        date < '2018-04-01' ~ '2',
                                        date < '2018-10-01' ~ '3',
                                        date < '2019-04-01' ~ '4',
                                        date < '2019-10-01' ~ '5',
                                        date < '2020-04-01' ~ '6')) %>%
      group_by(time_value, time_value_num, arm) %>%
      summarize(person_months = sum(person_months),
                n_cases = sum(n_cases)) %>%
      mutate(incidence_per_1000pm = n_cases / person_months * 1000,
             time_value = as.character(time_value),
             time_unit = agg_unit)
    
    inci_overall <- monthly_incidence_trial %>%
      group_by(arm) %>%
      summarize(person_months = sum(person_months),
                n_cases = sum(n_cases)) %>%
      mutate(incidence_per_1000pm = n_cases / person_months * 1000,
             time_value = 'Overall',
             time_value_num = 'Overall',
             time_unit = agg_unit)
    
    inci <- rbind(inci_annual, inci_overall) %>%
      mutate(time_value = factor(time_value, levels = c('June 2017-Sep 2017',
                                                        'Oct 2017-March 2018',
                                                        'April 2018-Sep 2018',
                                                        'Oct 2018-March 2019',
                                                        'April 2019-Sep 2019',
                                                        'Oct 2019-March 2020',
                                                        'Overall')))
  } else if (agg_unit == 'yearmonth'){
    inci <- monthly_incidence_trial %>%
      filter(date > '2017-05-01') %>%
      mutate(time_value = yearmonth, 
             time_value_num = yearmonth,
             time_unit = agg_unit)
  }
  
  # Pivot wider 
  inci_wide <- inci %>%
    dplyr::select(arm, time_value, time_value_num, time_unit, 
                  person_months, incidence_per_1000pm) %>%
    pivot_wider(
      names_from = arm,
      values_from = c(person_months, incidence_per_1000pm),
      id_cols = c(time_value, time_value_num, time_unit)
    )
  
  # Get IRRs
  inci_wide <- inci_wide %>%
    # mutate(incidence_per_1000pm_none = 40) %>%
    # Join the model predicted incidence for the no-intervention group to the trial data 
    left_join(model_inci_summary_wide_halfyear %>% dplyr::select(time_value, time_value_num, time_unit, incidence_none_median) %>%
                rename(incidence_per_1000pm_none = incidence_none_median)) %>%
    mutate(rtss_none_irr = (incidence_per_1000pm_rtss / incidence_per_1000pm_none),
           smc_none_irr = (incidence_per_1000pm_smc / incidence_per_1000pm_none),
           both_none_irr = (incidence_per_1000pm_both / incidence_per_1000pm_none),
           rtss_smc_irr = (incidence_per_1000pm_rtss / incidence_per_1000pm_smc),
           both_smc_irr = (incidence_per_1000pm_both / incidence_per_1000pm_smc),
           both_rtss_irr = (incidence_per_1000pm_both / incidence_per_1000pm_rtss),
           smc_rtss_irr = (incidence_per_1000pm_smc / incidence_per_1000pm_rtss)  )%>%
    mutate(expected_efficacy = 1 - (rtss_none_irr * smc_none_irr),
           ratio_pred_exp = (1-both_none_irr) / expected_efficacy,
           ratio_inci_rate_only = both_none_irr / (rtss_none_irr * smc_none_irr),
           
           inci_averted_model = incidence_per_1000pm_none - incidence_per_1000pm_both,
           cases_averted_model = inci_averted_model * 1000, # if pop is 1000
           
           inci_averted_expected = incidence_per_1000pm_none - 
             (incidence_per_1000pm_rtss * incidence_per_1000pm_smc)/incidence_per_1000pm_none,
           cases_averted_expected = inci_averted_expected * 1000, # if pop is 1000
           
           difference_inci_averted_pred_exp = inci_averted_model - inci_averted_expected,
           difference_cases_averted_pred_exp = cases_averted_model - cases_averted_expected)
  
  # Calculate efficacy
  inci_summary <- inci_wide %>%
    mutate(
      both_smc = 1- both_smc_irr,
      rtss_none =1-rtss_none_irr ,
      both_rtss = 1-both_rtss_irr,
      smc_none = 1-smc_none_irr  ,
      both_none = 1-both_none_irr,
      rtss_smc = 1-rtss_smc_irr  ,
      smc_rtss =1-smc_rtss_irr   )
  
  
  # Make long 
  inci_long <- inci_summary %>%
    dplyr::select(time_value, time_value_num, time_unit, starts_with('incidence')) %>%
    pivot_longer(cols = starts_with('incidence'),
                 names_to = c("arm", "stat"),
                 names_pattern = "incidence_(.*)_(.*)",
                 values_to = "value")%>%
    pivot_wider(
      names_from = stat,
      values_from = value
    )%>% 
    mutate(metric = 'incidence', comparison = NA)
  
  irrs_long <- inci_summary %>%
    dplyr::select(time_value, time_value_num, time_unit, both_smc:smc_rtss, expected_efficacy) %>%
    pivot_longer(cols =  c(both_smc:smc_rtss,expected_efficacy:expected_efficacy),
                 names_to = 'comparison',
                 values_to = "irr")%>%
    mutate(
      # Clean up the comparison names for better labels
      comparison = gsub("_", " vs ", comparison),
      comparison = ifelse(comparison == 'expected vs efficacy', 'Expected both vs none', comparison)
    ) %>%
    mutate(metric = 'efficacy',
           arm = NA)
  
  exp_ratio_long <- inci_summary %>%
    dplyr::select(time_value,time_value_num,  time_unit, contains('ratio'), 
                  -contains('averted'), -contains('ratio_inci_rate_only')) %>%
    pivot_longer(cols =  contains('ratio'),
                 names_to = 'statistic',
                 values_to = "ratio_pred_exp",
                 names_prefix = "ratio_pred_exp_") %>%
    pivot_wider(
      names_from = statistic,
      values_from = ratio_pred_exp
    ) %>% 
    mutate(metric = 'ratio pred to exp', arm = NA, comparison = NA)
  
  
  inci_averted_long <- inci_summary %>%
    dplyr::select(time_value,time_value_num,  time_unit, contains('averted'), contains('ratio_inci_rate_only')) %>%
    pivot_longer(cols =  c(contains('ratio_inci_rate_only'), contains('averted')),
                 names_to = c("metric"),
                 names_pattern = "^(.+)") %>%
    mutate(arm = NA, comparison = NA)
  
  inci_summary_all <- rbind(inci_long, 
                            irrs_long,
                            exp_ratio_long,
                            inci_averted_long)
  
  
  inci_summary <- inci_summary_all %>%
    mutate(shape_var = case_when(
      metric == 'efficacy' & grepl('Expected', comparison) ~ 'Expected',
      metric == 'efficacy' ~ 'Model-predicted',
      TRUE ~ NA)) %>%
    mutate(comparison = factor(comparison, levels = c("Expected both vs none", 'both vs none',
                                                      'both vs rtss','both vs smc',
                                                      'rtss vs none','smc vs none',
                                                      'rtss vs smc','smc vs rtss')))
  
  return(inci_summary)
}

#armfit = 3
trial_halfyearmali <- calc_irr_trial(agg_unit = 'halfyear',
                                     armfit = 3,
                                 country = 'Mali') %>%
  mutate(country = 'Mali')
trial_halfyearbf <- calc_irr_trial(agg_unit = 'halfyear',
                                   armfit = 3,
                                 country = 'BF') %>%
  mutate(country = 'BF')
trial_halfyear <- bind_rows(trial_halfyearmali, trial_halfyearbf)

inci_summary <- trial_halfyear

# Plot of the ratio of model-predicted to expected by aggregation unit 
pp <- ggplot(trial_halfyear %>% filter(metric == 'difference_inci_averted_pred_exp') %>%
               mutate(ypos = ifelse(value>= -0.2, value + 0.5, value-0.5)), 
       aes(x = time_value, y = value, 
           # fill = time_value, 
           group = time_value)) +
  # model estimated
  geom_col(fill = '#E15554') +
  geom_text(aes(x = time_value, y = ypos, label = round(value,2), group = value),
            size.unit = 'pt', size = 10) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
  labs(
    x = NULL,#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Half-year',
    y = "Difference in model-predicted trial observed versus expected\ncases averted per 1000 people of\ncombination vs no intervention",
    shape = NULL, linetype = NULL,
    color = NULL#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Study half-year'
  ) +
  theme_bw(base_size = 14) +
  facet_wrap(~country) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')

ggsave(filename = 'expected_trial_difference3.pdf', plot = pp, height = 6, width = 10)

# armfit = 2
trial_halfyearmali <- calc_irr_trial(agg_unit = 'halfyear',
                                     armfit = 2,
                                     country = 'Mali') %>%
  mutate(country = 'Mali')
trial_halfyearbf <- calc_irr_trial(agg_unit = 'halfyear',
                                   armfit = 2,
                                   country = 'BF') %>%
  mutate(country = 'BF')
trial_halfyear <- bind_rows(trial_halfyearmali, trial_halfyearbf)

inci_summary <- trial_halfyear

# Plot of the ratio of model-predicted to expected by aggregation unit 
pp <- ggplot(trial_halfyear %>% filter(metric == 'difference_inci_averted_pred_exp') %>%
               mutate(ypos = ifelse(value>= -0.2, value + 0.75, value-0.75)), 
             aes(x = time_value, y = value, 
                 # fill = time_value, 
                 group = time_value)) +
  # model estimated
  geom_col(fill = 'darkred' ) +
  geom_text(aes(x = time_value, y = ypos, label = round(value,2), group = value),
            size.unit = 'pt', size = 10) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'black') +
  labs(
    x = NULL,#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Half-year',
    y = "Difference in model-predicted trial observed versus expected\ncases averted per 1000 people of\ncombination vs no intervention",
    shape = NULL, linetype = NULL,
    color = NULL#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Study half-year'
  ) +
  theme_bw(base_size = 14) +
  facet_wrap(~country) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')
# pp
ggsave(filename = 'expected_trial_difference2.pdf', plot = pp, height = 6, width = 10)

# plotcompare <- ggplot(inci_summary %>% filter(metric == 'efficacy' & time_value == 'Overall'),
#        aes(x = comparison, y = irr,
#            shape = shape_var, fill = shape_var)) +
#   # model estimated
#   geom_col(position = position_dodge(width = 0.5), size = 1) +
#   scale_shape_manual(values = c('Expected' = 7, 'Model-predicted' = 16)) +
#   scale_color_manual(values = c('Expected' = "#8C6BB1", 'Model-predicted' ="#810F7C")) +
#   geom_point(data = model_inci_summary_all_halfyear %>% filter(metric == 'efficacy' & time_value == 'Overall'),
#              aes(x = comparison, y = median, fill = shape_var)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
#   labs(
#     x = "Arm Comparison",
#     y = "Relative efficacy (1-IRR)",
#     shape = NULL, linetype = NULL, fill = NULL,
#     color = NULL,#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Study half-year'
#     caption = 'Using model-predicted no-int incidence as comparison.\nBlack dots are from model output with fits from 23 Mar.\nBars are trial values.'#inci_wide$incidence_per_1000pm_none[1]
#   ) +
#   scale_y_continuous(breaks = seq(-1,1,0.2)) +
#   theme_minimal(base_size = 14) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# when 'none' incidence is 60, then expected and model predicted are the same
