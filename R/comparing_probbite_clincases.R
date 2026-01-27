# Compare the probability of a bite relative to the epidemic curve of clinical cases 
library(zoo)
library(survival)
library(survminer)
library(broom)
library(tidyverse)

source("R:/Kelly/synergy_orderly/shared/format_model_output.R")
source("R:/Kelly/synergy_orderly/shared/get_incidence.R")
source("R:/Kelly/synergy_orderly/shared/analyse_model_output.R")
source("R:/Kelly/synergy_orderly/shared/get_cox_efficacy.R")
source("R:/Kelly/synergy_orderly/shared/likelihood.R")

path <- 'R:/Kelly/synergy_orderly/src/sim_cohort_generic/outputs/'
foldername <- 'outputs_2025-12-01_2'
# sim_results <- readRDS(paste0(path,"outputs_2025-12-01_2/sim_results.rds"))
metadata_df <- readRDS(paste0(path, '/', foldername, "/metadata_df.rds"))

infectionrecords <-  readRDS(paste0(path,'/', foldername, "/infection_records.rds"))#purrr::map_df(sim_results, "infection_records")
params <-  readRDS(paste0(path,'/', foldername, "/parameter_df.rds"))#purrr::map_df(sim_results, 'params')
# parasitemia <- purrr::map_df(sim_results, 'parasitemia')
# infectionrecords <- filtered_df
# Format for incidence calculation 
formattedinfrecords <- lapply(params$sim_id, function(x){
  format_model_output(model_data = infectionrecords,
                      cohort = 'generic',
                      simulation = x)})
all <- bind_rows(formattedinfrecords)

# Calculate incidence 
testinci <- lapply(params$sim_id, function(x){
  aa <- all %>% filter(sim_id == x)

  get_incidence(df_children = metadata_df,
                          casedata = aa) %>%
    mutate(sim_id = x)
})
inci <- bind_rows(testinci)

# Summarize over all simulations 
inci_summ <- inci %>%
  group_by(arm, year, month, yearmonth) %>%
  summarize(across(c(incidence_per_1000pm),
                   list(lower = ~quantile(.x, 0.025, na.rm = TRUE),
                        median = ~quantile(.x, 0.5, na.rm = TRUE),
                        upper = ~quantile(.x, 0.975, na.rm = TRUE)),
                   .names = "{.col}_{.fn}") ) %>%
  # rename those variables with _median to be just the variable name 
  rename_with(.fn = \(x)sub("_median","", x)) 

# Plot incidence with p_bite overlaid 
pbite <- params$p_bite[[1]]
ggplot(inci_summ)+# %>% filter(sim_id == 'parameter_set_1_generic_FALSE')) +
  geom_line(aes(x = as.Date(yearmonth), y = incidence_per_1000pm, color = arm), linewidth = 2) +
  geom_point(aes(x = as.Date(yearmonth), y = incidence_per_1000pm, color = arm), size = 4) +
  geom_errorbar(aes(x = as.Date(yearmonth), ymin = incidence_per_1000pm_lower, ymax = incidence_per_1000pm_upper, color = arm),
                linewidth = 1, width = 20) +
  geom_point(data = pbite, aes(x = date, y = prob_lagged*1000), alpha = 0.2) +
  geom_vline(xintercept = as.Date(unlist(all$smc_dose_days[1][1:4]), origin = '2017-04-01'), linetype = 2, color = 'skyblue') +
  # facet_wrap(~sim_id) + 
  theme_bw() + 
  facet_wrap(~arm)

# analysis 
output <- all %>%
    filter(sim_id == 'parameter_set_1_generic_0.9')
  
model_output <- output %>%
  mutate(arm_smcref = factor(arm, levels = c('smc','rtss','both','none')),
         arm_rtssref = factor(arm, levels = c('rtss','smc','both','none')),
         arm_noneref = factor(arm, levels = c('none','rtss','smc','both')),
         country = 'BF',
         end_time = ifelse(start_time==end_time, end_time + 0.00001, end_time)) %>%
  filter(start_time != end_time) %>% ungroup()
  
kmsurvobj <- survival::survfit(Surv(start_time,
                                    end_time,
                                    event) ~ arm,#+ strata(country),
                               data = model_output)

cum_inci <- ggsurvplot(kmsurvobj, 
                       tables.theme = theme_cleantable(),
                       conf.int = TRUE,
                       fun = 'cumhaz',
                       ggtheme = theme_bw(base_size = 16),
                       palette = 'Dark2',
                       censor.size = 3)
cum_inci
survival <- ggsurvplot(kmsurvobj, 
                       tables.theme = theme_cleantable(),
                       conf.int = TRUE,
                       ggtheme = theme_bw(base_size = 16),
                       palette = 'Dark2',
                       censor.size = 3)
survival
# Infections per year 
infectionsperyear <- infectionrecords %>% filter(detectable==1 & 
                                                   sim_id == 'parameter_set_1_generic_0.9' &
                                                   detection_day < 365*3,
                                                 arm !='none') %>%
  group_by(rid) %>%
  summarize(n = n() / 3) %>%
  arrange(desc(n))%>%
  mutate(rid = factor(rid, levels = rid))
ggplot(infectionsperyear) + 
  geom_bar(aes(x = as.factor(rid), y = n), stat = 'identity', width = 1) + 
  theme_classic() + 
  labs(y = 'Average infections per year per person') + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
mean(infectionsperyear$n) # 7.2 --> with the EIR of 12, it is now 2.4 --> EIR 50 with treatment 3.3
median(infectionsperyear$n) # 4.3 --> with EIR of 12, 1.7 --> EIR 50 with treatment 2.7
# in the chandramohan trial it was mean of 1.8 and median of 1
infectionrecords %>% filter(detectable==1 & 
                              sim_id == 'parameter_set_1_generic_0.9' & 
                              detection_day < 365*3) %>%
  group_by(rid, arm) %>%
  summarize(n = n() / 3) %>%
  arrange(desc(n))%>%
  mutate(rid = factor(rid, levels = rid)) %>% group_by(arm) %>% summarize(median = median(n), mean = mean(n))

# Bites per year
bitesperyear <- infectionrecords %>% filter(sim_id == 'parameter_set_1_generic_0.9' & 
                              BSinfection_day < 365*3) %>%
  group_by(rid) %>%
  summarize(n = n() / 3) %>%
  arrange(desc(n)) %>%
  mutate(rid = factor(rid, levels = rid))
ggplot(bitesperyear) + 
  geom_bar(aes(x = as.factor(rid), y = n), stat = 'identity', width = 1) + 
  theme_classic() + 
  labs(y = 'Average infectious bites per year per person',
       caption = '*This is with treatment, which, in my model, removes people\nfrom the list of bitten children, so those are not counted here.') + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
mean(bitesperyear$n) #24.4 --> 7.8 with EIR 12 --> 12.5 as below
median(bitesperyear$n) # 17.7 --> 5.7 with EIR 12 --> 12.7 with EIR 50 but some removed because not susceptible

# Get efficacy
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
                      nonerefresults) %>%
  filter(term %in% c("Both vs. RTSS",'RTSS vs. SMC','Both vs. SMC'))

efficacy_plot <- ggplot(tidy_results)+
  geom_point(aes(x = term, y = VE)) +
  geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper), width = 0.2) +
  geom_hline(aes(yintercept = 0)) +
  labs(y = "Vaccine efficacy",
       x = NULL) +
  scale_y_continuous(breaks = seq(min(floor(tidy_results$VE_lower * 10)/10), 1.1, 0.2),
                     limits = c(min(floor(tidy_results$VE_lower * 10)/10),1)) +
  facet_wrap(~ year) + 
  theme_minimal()
efficacy_plot

# Monthly incidence
monthly_inci_trial <- readRDS("R:/Kelly/synergy_orderly/archive/trial_results/20250923-094333-2dd80ae3/monthly_incidence_trial.rds")
ll <- calculate_poisson_likelihood(trial_df = monthly_inci_trial, 
                                   model_df = inci %>% filter(sim_id == 'parameter_set_7_generic_0.9'))
ll
# Plot incidence for only 1 simulation 
incidence_plot <- inci %>% filter(sim_id == 'parameter_set_1_generic_0.9')  %>%
  # filter(intervention != 'none') %>%
  ggplot(aes(x = month, y = incidence_per_1000pm/10, color = arm)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower_per_1000/10, ymax = upper_per_1000/10, color = arm),
                alpha = 0.9, width = 10, linewidth = 1) +
  geom_line(linewidth = 1) +
  geom_line(data = monthly_inci_trial, aes(x = as.Date(yearmonth), y = incidence_per_1000pm), color = 'black', linetype = 2) +
  # geom_hline(aes(xintercept = c())) +
  scale_x_date(breaks = '3 months',
               labels = scales::label_date_short()) +
  labs(
    # title = "Monthly malaria incidence per 1000 person-months",
    x = "Year",
    y = "Incidence (per 1000 person-months)",
    color = "Study arm",
    fill = "Study arm"
  ) +
  theme_minimal(base_size = 16) + 
  facet_wrap(~arm)
incidence_plot


# Removal of infections that were within the last 7 days 
filtered_df <- infectionrecords %>%
  arrange(rid, detection_day) %>%
  group_by(rid) %>%
  mutate(prev_detection_day = lag(detection_day)) %>%
  filter(is.na(prev_detection_day) | (detection_day - prev_detection_day >= 7)) %>%
  select(-prev_detection_day) %>%
  ungroup()








# # look at infection records outside of fancy incidence calculation 
inc <- infectionrecords %>%
  mutate(detection_day = t - 50 + (threshold_day/2)) %>%
  group_by(detection_day, arm, sim_id) %>%
  # group_by(BSinfection_day, arm, sim_id) %>%
  summarize(ncases= sum(detectable)) %>%
  mutate(datedetection = as.Date(detection_day, origin = '2017-04-01'),
         # dateBSinfection = as.Date(BSinfection_day, origin = '2017-04-01'),
         month = lubridate::floor_date(datedetection, unit = 'month')) %>%
  group_by(month, arm, sim_id) %>%
  summarize(ncases = sum(ncases))

ggplot(inc) +
  geom_line(aes(x = month,
                y = ncases, color = sim_id)) +
  # geom_point(data = pbite, aes(x = date, y = prob_lagged*1000)) +
  geom_vline(xintercept = as.Date(unlist(all$smc_dose_days[1][1:4]), origin = '2017-01-01'), linetype = 1, color = 'grey50') +
  facet_wrap(~arm)

