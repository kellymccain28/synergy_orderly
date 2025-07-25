# Replicate survival analysis/calculation of incidence by month of trial results to model outputs
library(grates) # ifw e want to use greates_isoweek format like incidence does
library(survival)
library(tidyverse)
library(broom)
library(survminer)

orderly_strict_mode()

orderly_dependency(name = "sim_cohort", "latest()",
                   files = c("outputs/infection_records.rds", 
                             "outputs/parasitemia.rds", 
                             "outputs/metadata_children.rds",
                             "inputs.rds"))

orderly_resource(files = c('format_model_output.R'))
source('format_model_output.R')

orderly_shared_resource('get_cox_efficacy.R')
source('get_cox_efficacy.R')

dir.create("outputs/plots/")

orderly_artefact(description = 'plots and formatted model and analysis output',
                 files = c(
                   'outputs/plots/efficacy_plot.png',
                   'outputs/model_output_formatted.rds',
                   'outputs/plots/monthly_incidence_model.png',
                   'outputs/surv_analysis_model.rds', 
                   'outputs/monthly_incidence_model.rds'
                 ))

infection_records <- readRDS("outputs/infection_records.rds")
parasitemia <- readRDS("outputs/parasitemia.rds")
metadata_child <- readRDS("outputs/metadata_children.rds")
inputs <- readRDS("inputs.rds") # this is just to carry over the inputs 


# Format model output 
output <- format_model_output(infection_records)
saveRDS(output, 'outputs/model_output_formatted.rds')


# Do the same analysis as we do on trial data in trial_results.R
model_output <- output %>%
  mutate(arm_smcref = factor(intervention, levels = c('smc','vax','vaxsmc','none')),
         arm_rtssref = factor(intervention, levels = c('vax','smc','vaxsmc','none')),
         arm_noneref = factor(intervention, levels = c('none','vax','smc','vaxsmc')),
         country = 'BF') %>%
  filter(start_time != end_time)


# Survival analysis to reproduce results from trial
## SMC comparator by year and overall ---- 
smcrefresults <- get_cox_efficacy(df = model_output, 
                                  ref = 'arm_smcref',
                                  model = TRUE)

## RTSS comparator by year and overall ---- 
rtssrefresults <- get_cox_efficacy(df = model_output, 
                                   ref = 'arm_rtssref',
                                   model = TRUE)

#None comparator by year and overall ---- 
nonerefresults <- get_cox_efficacy(df = model_output, 
                                   ref = 'arm_noneref',
                                   model = TRUE)

# Plot the vaccine efficacies 
tidy_results <- rbind(smcrefresults, 
                      rtssrefresults, 
                      nonerefresults)

saveRDS(tidy_results, 'outputs/surv_analysis_model.rds')


# Plot results from survival analysis 
survresults <- tidy_results_model %>%
  filter(term %in% c("Both vs. RTSS",'RTSS vs. SMC','Both vs. SMC'))

efficacy_plot <- ggplot(survresults)+
  geom_point(aes(x = term, y = VE)) + 
  geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper), width = 0.2) +
  geom_hline(aes(yintercept = 0)) + 
  labs(y = "Vaccine efficacy",
       x = NULL) +
  scale_y_continuous(breaks = seq(min(floor(survresults$VE_lower * 10)/10), 1, 0.1),
                     limits = c(min(floor(survresults$VE_lower * 10)/10),1))

ggsave(filename = 'outputs/plots/efficacy_plot.png', efficacy_plot, height = 6, width = 6)



# Find the monthly incidence ---- 

# First need to calculate the total person-time per month by intervention
# function below is the same as the one in trial_results except with different argument names 
make_child_months <- function(child_id, intervention, start_date, end_date) {
  month_seq <- seq(floor_date(start_date, "month"), floor_date(end_date, "month"), by = "month")
  tibble(
    child_id = child_id,
    intervention = intervention,
    month = floor_date(month_seq, "month")
  )
}

# expand all children 
person_months_df_model <- metadata_child %>%
  mutate(fu_end_date = ymd('2020-03-31'),
         start_date = ymd('2017-04-01')) %>%
  select(child_id, intervention, start_date, fu_end_date) %>%
  pmap_dfr( ~ make_child_months(..1, ..2, ..3, ..4))

# get back start and end date 
person_months_df2_model <- person_months_df_model %>%
  mutate(start_date = ymd('2017-04-01'),
         end_date = ymd('2020-03-31'))

# Calculate actual person-time per child-month
# for each month, we need to find max an dmin -- so if the start date is=or > the floor date of the month and < last day of that month, 
# then we take either the first day of the month or the start date. if the start date is > last da of taht month, then we take the 
# last day of that month (I think????)
# for the end date, i think it's fine 
person_months_df3_model <- person_months_df2_model %>%
  mutate(
    month_start = #if_else(start_date >= floor_date(month, "month") & start_date < ceiling_date(month, "month"),
      pmax(floor_date(month, "month"), start_date),
    # ceiling_date(month, "month")),
    month_end = pmin(ceiling_date(month, "month"), end_date),
    days_at_risk = time_length(interval(month_start, month_end), unit = 'days'),#as.numeric(difftime(month_end, month_start, units = "days")) + 1,
    person_months = time_length(interval(month_start, month_end), unit = 'months'),#days_at_risk / 30.44,  # approx days per month
    person_years = time_length(interval(month_start, month_end), unit = 'years'),#days_at_risk / 365.25
  )

# Summarize by calendar month and arm
persontime_bymonth_model <- person_months_df3_model %>%
  mutate(year = lubridate::year(month), 
         month_num = month(month)) %>%
  group_by(year, month_num, intervention) %>%
  summarise(
    person_months = sum(person_months),
    person_years = sum(person_years),
    .groups = "drop"
  )

# Calculate incidence per month ----
cases_by_month_model <- model_output %>%
  filter(poutcome == 1) %>%
  mutate(month = lubridate::month(detection_day),
         year = lubridate::year(end_date)) %>%
  group_by(intervention, month, year) %>%
  summarize(n_cases = sum(poutcome), .groups = 'drop') %>%
  rename(month_num = month) 

# join t person time data 
monthly_inci_model <- persontime_bymonth_model %>%
  left_join(cases_by_month_model, by = c('year','month_num','intervention')) %>%
  mutate(
    n_cases = replace_na(n_cases, 0),
    incidence_per_1000pm = (n_cases / person_months) * 1000,
    rate = n_cases / person_months, 
    se = sqrt(n_cases) / person_months,
    lower = (qchisq(0.025, 2 * n_cases) / 2) / person_months, 
    upper = (qchisq(0.975, 2 * n_cases + 1) / 2) / person_months,
    incidence_per_1000pm = rate * 1000,
    lower_per_1000 = lower * 1000,
    upper_per_1000 = upper * 1000,
  ) %>%
  mutate(monthyear = paste0(month_num,'-',year))%>%
  mutate(date = make_date(year, month_num, 1))

saveRDS(monthly_inci_model, 'outputs/monthly_incidence_model.rds')

incidence_plot <- monthly_inci_model  %>%
  # filter(intervention != 'none') %>%
  ggplot(aes(x = date, y = incidence_per_1000pm, color = intervention)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower_per_1000, ymax = upper_per_1000, color = intervention),
                alpha = 0.9, width = 10, linewidth = 1) +
  geom_line(linewidth = 1) + 
  # geom_hline(aes(xintercept = c())) +
  scale_x_date(breaks = '3 months',
               labels = scales::label_date_short()) + 
  # facet_wrap(~intervention, nrow = 4) +
  # scale_y_continuous(breaks = seq(0,150,25)) +
  labs(
    title = "Monthly malaria incidence per 1000 person-months",
    x = "Year",
    y = "Incidence (per 1000 person-months)",
    color = "Study intervention",
    fill = "Study intervention"
  ) +
  theme_minimal(base_size = 16)


ggsave(filename = 'outputs/plots/monthly_incidence_model.png', bg = 'white', incidence_plot, height = 6, width = 12)


# Number of infections per day # incidence
ints <- c('smc','vax','vaxsmc','none')
weeks <- seq.Date(as.Date('2017-04-02'), as.Date('2020-03-31'), 7)
allweeks <- expand.grid(intervention = ints, 
                        week = weeks)# %>%
  # mutate(week = rep(1:length(weeks), each = length(ints)))

weeklyinf <- infection_records %>%
  mutate(week = floor_date(as.Date(detection_day, origin = '2017-04-01'), unit = 'week')) %>%
  merge( allweeks, all.y = TRUE) %>%
  group_by(intervention, week) %>%
  count()

n_infections <- ggplot(weeklyinf) + 
  geom_line(aes(x = week, y = n, color = intervention)) +
  geom_line(aes(x = week, y = n, color = intervention)) +
  facet_wrap(~ intervention) +
  labs(y = 'Weekly number of infections',
       x = 'Week since start of follow up period',
       caption = 'Assuming that the liver stage lasts 8 days') +
  theme(legend.position = 'none') + 
  theme_bw(base_size = 15)

ggsave(filename = 'outputs/plots/n_infections_model.png', bg = 'white', n_infections, height = 6, width = 12)


ggplot(parasitemia %>% filter(day1_BSinfection == 19 & intervention == 'smc')) + 
  geom_line(aes(x=time_ext/7, y = prob_smckill))

# Proportion detectable ---- 
prop_det <- ggplot(infection_records ) + 
  geom_bar(aes(x = intervention, group = as.factor(detectable), fill = as.factor(detectable)),
           position = 'fill') + 
  scale_fill_manual(values = c('darkmagenta','goldenrod')) +
  labs(#title = 'Detectable infections in each intervention group',
    fill = 'Detectable')

ggsave(filename = 'outputs/plots/prop_detectable_model.png', prop_det, height = 6, width = 6)



# Cumulative incidence ---- 
output <- output %>%
  mutate(end_time = ifelse(start_time==end_time, end_time + 0.0001, end_time),
         intervention = factor(intervention, levels = c('vaxsmc','vax','smc','none')))

kmsurvobj <- survfit(Surv(start_time, 
                          end_time, 
                          event) ~ intervention,#+ strata(country), 
                     data = output)

cum_inci <- ggsurvplot(kmsurvobj, #group.by = "country",
           # linetype = "strata",
           # facet.by = 'country',
           tables.theme = theme_cleantable(),
           conf.int = TRUE,
           fun = 'cumhaz',
           # risk.table = TRUE,
           # cumevents = TRUE,
           ggtheme = theme_bw(base_size = 16),
           palette = 'Dark2',
           censor.size = 3)

ggsave(filename = 'outputs/plots/cum_inci_model.png', cum_inci$plot, height = 6, width = 10)


survival <- ggsurvplot(kmsurvobj, #group.by = "country",
                       # linetype = "strata",
                       # facet.by = 'country',
                       tables.theme = theme_cleantable(),
                       conf.int = TRUE,
                       # fun = 'cumhaz',
                       # risk.table = TRUE,
                       # cumevents = TRUE,
                       ggtheme = theme_bw(base_size = 16),
                       palette = 'Dark2',
                       censor.size = 3)

ggsave(filename = 'outputs/plots/survival_model.png', survival$plot, height = 6, width = 10)


