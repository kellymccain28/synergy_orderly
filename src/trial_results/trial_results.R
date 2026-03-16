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
                             'data/delivery_detail.rds',
                             'data/weekly.rds'))

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
ggsave(filename = 'efficacy_trial.pdf', efficacies, height = 5, width = 8)


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

ggplot(results )+#%>% filter(term != 'factor(arm)rtss'))+
  geom_point(aes(x = term, y = VE)) + 
  geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper), width = 0.2) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  labs(x = '') +
  theme_bw(base_size = 16) 

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
  select(date, arm) %>% 
  mutate(color = '#709176') 
smc_lines_Mali <-  smc_dates %>% ungroup() %>%
  filter(country == 'Mali' & arm != 'rtss') %>%
  select(date, arm) %>% 
  mutate(color = '#709176') 

vaxdates <- readRDS('R:/Kelly/synergy_orderly/shared/median_rtss_dates.rds')
rtss_lines_BF <-  vaxdates %>% ungroup() %>%
  filter(country == 'BF' & arm != 'smc') %>%
  select(date, arm) %>% 
  mutate(color = '#59114D') 
rtss_lines_Mali <-  vaxdates %>% ungroup() %>%
  filter(country == 'Mali' & arm != 'smc') %>%
  select(date, arm) %>% 
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
monthlyincidenceplot
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

# Average delivery times for each intervention / country 
ggplot(delivery %>% filter(arm !='rtss')) +
  geom_histogram(aes(x = nsmc_received, group = arm, fill = arm),
                 position = 'dodge') + 
  facet_wrap(~ country) +
  scale_x_continuous(breaks =seq(0,12)) + 
  scale_fill_manual(values = mycols) +
  labs(x = 'Number of SMC rounds received (maximum of 12)',
       y = 'Number of people',
       fill = 'Intervention arm') + 
  theme_bw(base_size =  14)
ggsave('nsmc_received_bycountryandarm.pdf', plot = last_plot(), height = 5, width = 10)

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
  group_by(country, arm) %>%
  summarize(across(contains('date'),
                   list(median = ~median(.x, na.rm = TRUE)),
                   .names = "{.col}_{.fn}")) #%>%
# pivot_longer(cols = v1_date_median:fu_end_date_median)
vax_dates_avg <- delivery_avg %>%
  dplyr::select(country, arm, v1_date_median, v2_date_median, v3_date_median, boost1_date_median, boost2_date_median) %>%
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
  facet_wrap(~country + arm, scales = 'free')

median(delivery$y1p1d1_date_received, na.rm = TRUE) # 7-27
median(delivery$y1p2d1_date_received, na.rm = TRUE) # 8-24 
median(delivery$y1p3d1_date_received, na.rm = TRUE) # 9-23

smcdates <- delivery_avg %>%
  dplyr::select(country, arm, contains('d3')) %>%
  pivot_longer(cols = contains('date_received'),
               names_to = 'smcdose',
               values_to = 'date') 
ggplot(smcdates) +
  geom_vline(aes(xintercept = date, color = smcdose, linetype = arm)) + facet_wrap(~country)

saveRDS(smcdates, 'R:/Kelly/synergy_orderly/shared/median_smc_dates.rds')


# Parasitaemia 
weekly <- weekly %>%
  mutate(poutcome= ifelse(pf_asex_fdensity >= 5000, '1', '0'))

ggplot(weekly) +
  geom_jitter(aes(x = arm, y = pf_asex_fdensity, color = poutcome, shape = country)) + 
  geom_violin(aes(x = arm, y= pf_asex_fdensity), alpha = 0.4) +
  scale_y_log10() +
  theme_classic()

ggplot(weekly) +
  geom_point(aes(x = dateweekly, y = pf_asex_fdensity, color = arm)) + 
  scale_y_log10() +
  theme_classic()

weekly %>% 
  group_by(arm) %>%
  summarise(meanpb= mean(pf_asex_fdensity, na.rm = TRUE),
            medianpb = median(pf_asex_fdensity, na.rm = TRUE))


# Calculate IRRs as in model outputs 
source("R:/Kelly/synergy_orderly/src/make_figures_modeldev/bootstrap_metric.R")
library(purrr)
model_inci_summary_all_halfyear <- readRDS("R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs/outputs_2026-03-02/inci_summary_all_halfyear.rds") %>%
  mutate(shape_var = case_when(
    metric == 'efficacy' & grepl('Expected', comparison) ~ 'Expected',
    metric == 'efficacy' ~ 'Model-predicted',
    TRUE ~ NA)) 
model_inci_summary_wide_halfyear <- readRDS("R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs/outputs_2026-03-02/inci_summary_wide_halfyear.rds")

inci <- monthly_inci# trial
agg_unit = 'halfyear'
inci_annual <- inci %>%
  filter(!is.na(date)) %>%
  filter(date > '2017-05-01') %>%
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

inci_overall <- inci %>%
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
  left_join(model_inci_summary_wide_halfyear %>% select(time_value, time_value_num, time_unit, incidence_none_median) %>%
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

# Calculate median and IQR for each metric by time aggregation unit 
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

ggplot(inci_summary %>% filter(metric == 'efficacy'), aes(x = comparison, y = irr, 
                                                          color = time_value, group = time_value, shape = shape_var)) +
  # model estimated
  geom_point(position = position_dodge(width = 0.5), size = 1) +
  scale_shape_manual(values = c('Expected' = 7, 'Model-predicted' = 16)) + 
  # scale_linetype_manual(values = c('Expected: both vs none' = 2)) +
  # scale_color_brewer(palette = 'BuPu') +
  # scale_color_manual(values = colors) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
  labs(
    x = "Arm Comparison",
    y = "Relative efficacy (1-IRR)",
    # title = "Median IRR with 95% CI by Intervention Comparison",
    shape = NULL, linetype = NULL,
    color = if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Study half-year'
  ) +
  scale_y_continuous(breaks = seq(-1,1,0.2)) + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(inci_summary %>% filter(metric == 'efficacy' & time_value == 'Overall'), 
       aes(x = comparison, y = irr, 
           shape = shape_var, fill = shape_var)) +
  # model estimated
  geom_col(position = position_dodge(width = 0.5), size = 1) +
  scale_shape_manual(values = c('Expected' = 7, 'Model-predicted' = 16)) + 
  scale_color_manual(values = c('Expected' = "#8C6BB1", 'Model-predicted' ="#810F7C")) + 
  geom_point(data = model_inci_summary_all_halfyear %>% filter(metric == 'efficacy' & time_value == 'Overall'), 
             aes(x = comparison, y = median, fill = shape_var)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
  labs(
    x = "Arm Comparison",
    y = "Relative efficacy (1-IRR)",
    shape = NULL, linetype = NULL, fill = NULL,
    color = NULL,#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Study half-year'
    caption = 'Model-predicted no int incidence\nBlack dots are from model output with fits from 2 Mar'#inci_wide$incidence_per_1000pm_none[1]
  ) +
  scale_y_continuous(breaks = seq(-1,1,0.2)) + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# when 'none' incidence is 60, then expected and model predicted are the same

# Plot of the ratio of model-predicted to expected by aggregation unit 
ggplot(inci_summary %>% filter(metric == 'difference_inci_averted_pred_exp'), 
       aes(x = time_value, y = value, 
           color = time_value, group = time_value)) +
  # model estimated
  geom_point(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
  labs(
    x = NULL,#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Half-year',
    y = "Difference in model-predicted versus expected\ncases averted per 1000 people of\ncombination vs no intervention",
    shape = NULL, linetype = NULL,
    color = NULL#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Study half-year'
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')
