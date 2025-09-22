# Reproducing trial outputs from Chandramohan et al 2021
library(tidyverse)
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
tidy_results <- rbind(smcrefresults, rtssrefresults)

saveRDS(tidy_results, 'surv_analysis_trial.rds')

# Get efficacy for 3 versus 2 doses 
df <- primary %>%
  # filter(nprimary %in% c(3,2)) %>%
  mutate(nprimary = factor(nprimary, levels = c(3, 2, 1))) %>%
  filter(arm == 'both'| arm == 'rtss') %>%#  & year == 1)
  filter(arm == 'rtss')
coxdoses <- coxph(Surv(start_time, end_time, event) ~ factor(nprimary),# + factor(country),
                  data = df,
                  cluster = rid,
                  ties = "efron")
results <- tidy(coxdoses,
                exponentiate = TRUE,
                conf.int = TRUE) %>%
  # Calculate vaccine efficacy and its confidence intervals
  mutate(VE = (1 - estimate) * 100,              # VE = 1 - HR
    VE_lower = (1 - conf.high) * 100,       # Lower CI for VE = 1 - Upper CI for HR
    VE_upper = (1 - conf.low) * 100        # Upper CI for VE = 1 - Lower CI for HR
  ) %>% mutate(
    n_events = coxdoses$nevent,
    n_obs = coxdoses$n
  )
ggplot(results)+
  geom_point(aes(x = term, y = VE)) + 
  geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper), width = 0.2) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  labs(x = '') +
  theme_bw(base_size = 16) #+ labs(caption = 'rtss only')

# Plot efficacy
efficacies <- ggplot(tidy_results %>% filter(term %in% c("Both vs. RTSS",'RTSS vs. SMC','Both vs. SMC')))+
  geom_point(aes(x = term, y = VE)) + 
  geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper), width = 0.2) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_y_continuous(breaks = seq(-0.4, 1, 0.1),
                     limits = c(-0.3, 1)) +
  facet_wrap(~year) + 
  labs(x = '') +
  theme_bw(base_size = 16)

ggsave(filename = 'efficacy_trial.png', efficacies, height = 8, width = 8)


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
  ggplot(aes(x = date, y = incidence_per_1000pm, color = arm)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower_per_1000, ymax = upper_per_1000, color = arm),
                alpha = 0.9, width = 10, linewidth = 1) +
  geom_line(aes(x = date, y = incidence_per_1000pm, color = arm), linewidth = 1) +
  # facet_wrap(~arm, nrow = 3) +
  scale_y_continuous(breaks = seq(0,150,25)) +
  scale_x_date(breaks = '3 months',
               labels = scales::label_date_short()) + 
  labs(
    title = "Monthly malaria incidence per 1000 person-months",
    x = "Month",
    y = "Incidence (per 1000 person-months)",
    color = "Study Arm",
    fill = "Study Arm"
  ) +
  theme_minimal(base_size = 14)

ggsave("trial_monthlyincidence.png", plot = monthlyincidenceplot, bg = 'white', width = 8, height = 6)

















 ###### Get follow-up time among cases and non-cases-----
# make_child_months <- function(rid, arm, v1_date, fu_end_date) {
#   month_seq <- seq(floor_date(v1_date, "month"), floor_date(fu_end_date, "month"), by = "month")
#   tibble(
#     rid = rid,
#     arm = arm,
#     month = floor_date(month_seq, "month")
#   )
# }
# 
# # expand all children 
# person_months_df <- children %>%
#   mutate(fu_end_date = ymd('2020-03-31')) %>%
#   select(rid, arm, v1_date, fu_end_date) %>%
#   pmap_dfr( ~ make_child_months(rid = ..1, arm = ..2, v1_date = ..3, fu_end_date = ..4))
# 
# # join back v1date and end date 
# person_months_df <- person_months_df %>%
#   left_join(children %>% mutate(fu_end_date = ymd('2020-03-31')) %>%
#               select(rid, v1_date, fu_end_date), by = 'rid')
# 
# # Calculate actual person-time per child-month
# person_months_df <- person_months_df %>%
#   mutate(
#     month_start = pmax(floor_date(month, "month"), v1_date),
#     month_end = pmin(ceiling_date(month, "month") - days(1), fu_end_date),
#     days_at_risk = time_length(interval(month_start, month_end), unit = 'days'),#as.numeric(difftime(month_end, month_start, units = "days")) + 1,
#     person_months = time_length(interval(month_start, month_end), unit = 'months'),#days_at_risk / 30.44,  # approx days per month
#     person_years = time_length(interval(month_start, month_end), unit = 'years'),#days_at_risk / 365.25
#   )
# 
# # Summarize by calendar month and arm
# persontime_bymonth <- person_months_df %>%
#   mutate(year = lubridate::year(month), 
#          month = month(month),
#          yearmonth = zoo::as.yearmon(month)) %>%
#   group_by(year, month,yearmonth, arm) %>%
#   summarise(
#     person_months = sum(person_months),
#     person_years = sum(person_years),
#     .groups = "drop"
#   )
# 
# # # Summarize by calendar year and arm 
# # persontime_byyear <- person_months_df %>%
# #   mutate(year = lubridate::year(month), 
# #          month_num = month(month)) %>%
# #   group_by(year, arm) %>%
# #   summarise(
# #     person_months = sum(person_months),
# #     person_years = sum(person_years),
# #     .groups = "drop"
# #   )
# # 
# # # Summarize by arm 
# # persontime_overall <- person_months_df %>%
# #   mutate(year = lubridate::year(month), 
# #          month_num = month(month)) %>%
# #   group_by(arm) %>%
# #   summarise(
# #     person_months = sum(person_months),
# #     person_years = sum(person_years),
# #     .groups = "drop"
# #   )
# 
# # per protocol population for each trial year should also be those that received all doses of vaccine and attneded all 4 chemoprevention visits in that year 
# # primary_pp <- primary %>%
# #   filter(!is.na(v1_date) & !is.na(v2_date) & !is.na(v3_date)) 
# 
# 
# # Calculate person time (years and months at risk)
# # persontime_bymonth <- mitt %>%
# #   group_by(month, year, arm) %>%
# #   mutate(person_months_per_person = time_length(interval(v1_date, ymd('2020-03-31')), unit = 'months'),
# #          person_years_per_person = time_length(interval(v1_date, ymd('2020-03-31')), unit = 'year')) %>%
# #   summarize(person_months = sum(person_months_per_person),
# #             person_years = sum(person_years_per_person), .groups = 'drop')
# # 
# # persontime <- mitt %>%
# #   ungroup() %>%
# #   group_by(arm) %>%
# #   mutate(person_months_per_person = time_length(interval(v1_date, ymd('2020-03-31')), unit = 'months'),
# #          person_years_per_person = time_length(interval(v1_date, ymd('2020-03-31')), unit = 'year')) %>%
# #   summarize(person_months = sum(person_months_per_person),
# #             person_years = sum(person_years_per_person))
# # 
# # events_by_month <- mitt %>%
# #   group_by(arm, month, year) %>%
# #   summarize(n_events = sum(poutcome), .groups = 'drop') %>%
# #   left_join(persontime_bymonth, by = c('arm','month','year')) %>%
# #   mutate(incidence = n_events / person_months * 1000) %>%
# #   arrange(year, month) %>%
# #   mutate(month_index = row_number(),
# #          date = factor(paste0(month, "-",year)))
# 
# # Calculate incidence per month ----
# cases_by_month <- mitt %>%
#   filter(poutcome == 1 )%>%
#   mutate(yearmonth = zoo::as.yearmon(month)) %>%
#   group_by(arm, month, year, yearmonth) %>%
#   summarize(n_cases = sum(poutcome), .groups = 'drop')
# 
# # join t person time data 
# monthly_inci <- persontime_bymonth %>%
#   left_join(cases_by_month, by = c('yearmonth','year','month','arm')) %>%
#   mutate(
#     n_cases = replace_na(n_cases, 0),
#     incidence_per_1000pm = (n_cases / person_months) * 1000,
#     rate = n_cases / person_months, 
#     se = sqrt(n_cases) / person_months,
#     lower = (qchisq(0.025, 2 * n_cases) / 2) / person_months, 
#     upper = (qchisq(0.975, 2 * n_cases + 1) / 2) / person_months,
#     incidence_per_1000pm = rate * 1000,
#     lower_per_1000 = lower * 1000,
#     upper_per_1000 = upper * 1000,
#   ) 
# 
# monthly_inci %>%
#   mutate(date = make_date(year, month, 1)) %>%
#   ggplot(aes(x = date, y = incidence_per_1000pm, color = arm)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = lower_per_1000, ymax = upper_per_1000, color = arm),
#                 alpha = 0.9, width = 10, linewidth = 1) +
#   geom_line(aes(x = date, y = incidence_per_1000pm, color = arm), linewidth = 1) +
#   # facet_wrap(~arm, nrow = 3) +
#   scale_y_continuous(breaks = seq(0,150,25)) +
#   scale_x_date(breaks = '3 months',
#                labels = scales::label_date_short()) + 
#   labs(
#     title = "Monthly malaria incidence per 1000 person-months",
#     x = "Month",
#     y = "Incidence (per 1000 person-months)",
#     color = "Study Arm",
#     fill = "Study Arm"
#   ) +
#   theme_minimal(base_size = 16)
# 
# 
# 
# # Calculate incidence per month in first infection df----
# cases_by_month <- mitt_first %>%
#   filter(poutcome == 1) %>%
#   group_by(arm, month, year) %>%
#   summarize(n_cases = sum(poutcome), .groups = 'drop') %>%
#   rename(month_num = month)
# 
# cases_by_month %>%
#   mutate(date = make_date(year, month_num, 1)) %>%
#   ggplot(aes(x = date, y = n_cases, color = arm)) +
#   geom_point() +
#   facet_wrap(~arm, nrow = 3) +
#   scale_y_continuous(breaks = seq(0,150,25)) +
#   labs(
#     title = "Monthly malaria incidence per 1000 person-months",
#     x = "Month",
#     y = "Incidence (per 1000 person-months)",
#     color = "Study Arm",
#     fill = "Study Arm"
#   ) +
#   theme_minimal()
# 
# # Kaplan-Meier model to estimate survival among first infections only
# kmsurvobj <- survfit(Surv(time_to_event, event) ~ arm, 
#                      data = mitt_first)
# plot(kmsurvobj, mark.time = TRUE, 
#      fun = 'cumhaz',
#      conf.int = TRUE, 
#      col = c('red','blue','purple'))
# survdiff(Surv(time_to_event, event) ~ arm, 
#          data = mitt_first)
# 
# 
# # a bit confused because the cumulative hazard curves look super different compared to the one published in the paper -- maybe because the FU times are incorrect?? 
# 
# # Parametric survival analysis
# # https://bookdown.org/drki_musa/dataanalysis/parametric-survival-analysis.html
# weib_model <- survreg(Surv(time_to_event, event) ~ arm,
#                       data = mitt_first, 
#                       dist = 'weibull')
# summary(weib_model)
# # this shows both the hazard ratio (HR) for prop hazards, and the event time ratio (AFT modl)
# ConvertWeibull(weib_model, conf.level = 0.95)
# 
# weib_model2 <- flexsurvreg(Surv(time_to_event, event) ~ arm, 
#                            data = mitt_first, 
#                            dist = 'weibull')
# weib_model2
# WeibullDiag(Surv(time = mitt_first$time_to_event, event = mitt_first$event == 1) ~ arm, 
#             data = mitt_first)
# # # Fit Poisson model (rates by arm and calendar month)
# # poisson_model <- glm(
# #   n_events ~ arm + date,  # or use month_index as numeric trend
# #   offset = log(person_months),
# #   family = poisson(link = "log"),
# #   data = events_by_month
# # )
# # summary(poisson_model)
# # 
# # # Predict on the original data
# # # pred_log <- predict(poisson_model, type = "link", se.fit = TRUE, offset = log(person_months))
# # pred_log <- augment(poisson_model,
# #                     type.predict = c('link'),
# #                     se_fit = TRUE)
# # 
# # # Convert back to incidence rate scale
# # events_by_month <- events_by_month %>%
# #   mutate(
# #     rate = exp(pred_log$.fitted) *1000,  # per 1000 person-months
# #     lower = exp(pred_log$.fitted - 1.96 * pred_log$.se.fit)  *1000,
# #     upper = exp(pred_log$.fitted + 1.96 * pred_log$.se.fit) *1000
# #   )
# # 
# # ggplot(events_by_month, aes(x = month_index, y = rate/1000, color = arm)) +
# #   # geom_line() +
# #   geom_point() +
# #   geom_errorbar(aes(ymin = lower/1000, ymax = upper/1000)) +
# #   facet_wrap(~ arm, nrow = 3) +
# #   ylab("Malaria Incidence (per 1000 person-months)") +
# #   xlab("Month") +
# #   theme_minimal()
# 
# 
# # other covariates
# cox <- coxph(Surv(start, stop, event) ~ interaction(arm,month) + strata(country), 
#              cluster = rid,
#              data = mitt_first)
# coxsumm <- summary(cox)
# resultTable <- cbind(coxsumm$coefficients, coxsumm$conf.int[, c('lower .95', 'upper .95')])
# round(resultTable, 4)
# # 1.18(1.07, 1.30) is the hazard ratio of clinical case of malaria in rtss arm versus smc arm, and 
# # 0.73(0.65, 0.82) is hazard ratio of clincal cases of malaria in both arm versus smc arm
# 
# #https://yzwisalaity.github.io/Baseline-Hazard-Plot/ or Nelson-Aalen
# bashaz <- basehaz(cox)
# ggplot(bashaz) + 
#   geom_line(aes(x = time, y = hazard, color = strata))
# 
# # Cox model 
# # to get hazards
# cox_model <- coxph(
#   Surv(start, stop, event) ~ arm + cluster(rid) + strata(country) ,
#   data = mitt %>% mutate(arm = factor(arm, levels =c('smc','rtss','both'))),
#   robust = TRUE
# )
# 
# # Get hazard ratios with 95% CIs
# summary(cox_model)$conf.int %>%
#   as.data.frame() %>%
#   select(HR = `exp(coef)`, Lower = `lower .95`, Upper = `upper .95`) %>%
#   mutate(
#     ProtectiveEfficacy = (1 - HR) * 100,
#     PE_Lower = (1 - Upper) * 100,
#     PE_Upper = (1 - Lower) * 100
#   )
# 
# # Get Nelson-Aalen cumulative ahzard estimate
# fit <- survfit(Surv(start, stop, event) ~ arm, 
#                data = mitt %>% mutate(arm = factor(arm, levels =c('smc','rtss','both'))))
# 
# plot(fit, fun = "cumhaz", col = 1:3, lty = 1:3,
#      xlab = "Time", ylab = "Cumulative hazard",
#      main = "Nelson-Aalen cumulative hazard by arm")
# legend("topleft", legend = levels(primary_pp$arm), col = 1:3, lty = 1:3)




####
