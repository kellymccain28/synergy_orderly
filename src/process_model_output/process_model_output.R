# Replicate survival analysis/calculation of incidence by month of trial results to model outputs
library(grates) # ifw e want to use greates_isoweek format like incidence does
library(survival)
library(tidyverse)
library(broom)

orderly_strict_mode()

orderly_dependency(name = "sim_cohort", "latest()",
                   files = c("outputs/infection_records.rds", 
                             "outputs/parasitemia.rds", 
                             "outputs/metadata_children.rds",
                             "inputs.rds"))

orderly_resource(files = c('format_model_output.R'))
source('format_model_output.R')

orderly_artefact(description = 'plots and formatted model and analysis output',
                 files = c(
                   'efficacy_plot.png',
                   'model_output_formatted.png',
                   'monthly_incidence_model.png',
                   'surv_analysis_model.png', 
                   'monthly_incidence_model.rds'
                 ))

infection_records <- readRDS("outputs/infection_records.rds")
parasitemia <- readRDS("outputs/parasitemia.rds")
metadata_child <- readRDS("outputs/metadata_children.rds")
inputs <- readRDS("inputs.rds") # this is just to carry over the inputs 


# Format model output 
output <- format_model_output(infection_records)
saveRDS(output, 'model_output_formatted.rds')


# Do the same analysis as we do on trial data in trial_results.R
model_output <- output %>%
  mutate(arm_smcref = factor(intervention, levels = c('smc','vax','vaxsmc','none')),
         arm_rtssref = factor(intervention, levels = c('vax','smc','vaxsmc','none')),
         country = 'BF') %>%
  filter(start_time != end_time)


# Survival analysis to reproduce results from trial
# with Smc as comparator :
ag_model_smcref_model <- coxph(Surv(start_time, 
                                    end_time, 
                                    event) ~ factor(arm_smcref),# + factor(country), 
                               data = model_output,
                               cluster = child_id,
                               ties = "efron")

# Display results
summary(ag_model_smcref_model)

coxsumm_model<-as.data.frame(summary(ag_model_smcref_model)$coefficients)

tidy_results_smcref_model <- tidy(ag_model_smcref_model,
                                  exponentiate = TRUE,
                                  conf.int = TRUE) %>%
  # Calculate vaccine efficacy and its confidence intervals
  mutate(
    term = case_when(
      term == 'factor(arm_smcref)vax' ~ "RTSS vs. SMC",
      term == 'factor(arm_smcref)vaxsmc' ~ "Both vs. SMC",
      term == 'factor(country)Mali' ~ 'Mali vs. BF',
      term == 'factor(arm_smcref)none' ~ "None vs. SMC",
      TRUE ~ term
    ),
    VE = 1 - estimate,              # VE = 1 - HR
    VE_lower = 1 - conf.high,       # Lower CI for VE = 1 - Upper CI for HR
    VE_upper = 1 - conf.low         # Upper CI for VE = 1 - Lower CI for HR
  )

# tidy_results_smcref_model

# with rtss as comparator: 
ag_model_rtssref_model <- coxph(Surv(start_time, 
                                     end_time, 
                                     event) ~ factor(arm_rtssref),# + factor(country), 
                                data = model_output,
                                cluster = child_id,
                                ties = "efron")

# Display results
# summary(ag_model_rtssref_model)

coxsumm2_model<-as.data.frame(summary(ag_model_rtssref_model)$coefficients)

tidy_results_rtssref_model <- tidy(ag_model_rtssref_model,
                                   exponentiate = TRUE,
                                   conf.int = TRUE) %>%
  # Calculate vaccine efficacy and its confidence intervals
  mutate(
    term = case_when(
      term == 'factor(arm_rtssref)smc' ~ "SMC vs. RTSS",
      term == 'factor(arm_rtssref)vaxsmc' ~ "Both vs. RTSS",
      term == 'factor(country)Mali' ~ 'Mali vs. BF',
      term == 'factor(arm_rtssref)none' ~ "None vs. RTSS",
      TRUE ~ term
    ),
    VE = 1 - estimate,              # VE = 1 - HR
    VE_lower = 1 - conf.high,       # Lower CI for VE = 1 - Upper CI for HR
    VE_upper = 1 - conf.low         # Upper CI for VE = 1 - Lower CI for HR
  ) 

# tidy_results_rtssref_model


# Plot the vaccine efficacies 
tidy_results_model <- rbind(tidy_results_rtssref_model, tidy_results_smcref_model)

saveRDS(tidy_results_model, 'surv_analysis_model.rds')

# Plot results from survival analysis 
efficacy_plot <- ggplot(tidy_results_model %>% filter(term %in% c("Both vs. RTSS",'RTSS vs. SMC','Both vs. SMC')))+
  geom_point(aes(x = term, y = VE)) + 
  geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper), width = 0.2) +
  geom_hline(aes(yintercept = 0)) + 
  labs(y = "Vaccine efficacy",
       x = NULL) +
  scale_y_continuous(breaks = seq(0, 1, 0.1),
                     limits = c(0,1))

ggsave(filename = 'efficacy_plot.png', efficacy_plot, height = 6, width = 6)



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

saveRDS(monthly_inci_model, 'monthly_incidence_model.rds')

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


ggsave(filename = 'monthly_incidence_model.png', incidence_plot, height = 6, width = 12)












