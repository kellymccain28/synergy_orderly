library(haven)
library(ggplot2)
library(dplyr)
library(janitor)
library(gtsummary)
library(wesanderson)
library(tidyr)
library(stringr)
library(survival)
library(survminer)
library(broom)
library(gtsummary)
library(cyphr)
library(orderly)
library(lubridate)

key <- cyphr::data_key()

orderly_strict_mode()

orderly_artefact(description = 'encrypted cleaned files',
                 files = c('data/primary_persontime.rds',
                           'data/mitt.rds',
                           'data/weekly.rds',
                           'data/serology.rds',
                           'data/delivery_detail.rds',
                           'data/primary.rds',
                           'data/children.rds'))

pathtodata = "C:/Users/kem22/OneDrive - Imperial College London/PhD Admin/Within host/RTSS SMC Data for Kelly_20250228/"

files <- list.files(paste0(pathtodata,'from Paul'),
                    full.names = TRUE)

data <- lapply(files, read_dta)

# Light data cleaning 
cleandata <- function(df){
  df <- df %>%
    mutate(sex = ifelse(sex == 1, 'm', 
                        ifelse(sex == 2, 'f', sex)),
           arm = case_when(
             arm == 1 ~ 'smc',
             arm == 2 ~ 'rtss',
             arm == 3 ~ 'both',
             TRUE ~ as.character(arm)
           ),
           country = case_when(
             grepl('B', rid) ~ 'BF',
             grepl('M', rid) ~ 'Mali',
             TRUE ~ NA
           ))
  
  return(df)
}

children_dob_raw <- as.data.frame(data[1]) %>%
  cleandata() 

children_raw <- as.data.frame(data[2])%>%
  cleandata()

primary_persontime <- as.data.frame(data[3])

primary_raw <- as.data.frame(data[4]) 

delivery_raw <- as.data.frame(data[5])%>%
  cleandata()%>%
  mutate(last_primary_vac = coalesce(v3_date, v2_date, v1_date))

delivery_detail_raw <- as.data.frame(data[6])%>%
  cleandata() %>%
  mutate(last_primary_vac = coalesce(v3_date, v2_date, v1_date))

serology_raw <- as.data.frame(data[7]) 

weekly1_raw <- as.data.frame(data[8]) 

weekly2_raw <- as.data.frame(data[9]) 

weekly3_raw <- as.data.frame(data[10]) 

children <- children_raw %>% 
  left_join(delivery_raw %>% 
              select(rid, arm, sex, 
                     v1_date, v2_date, v3_date, last_primary_vac, boost1_date, boost2_date, 
                     ends_with('date_received'), -contains('d2'), -contains('d3'), # d2 and d3 are the 2nd and third doses per month
                     nprimary, nsmc_received)) %>%
  mutate(dob = v1_date - age_months_v1 * 30) %>%
  # mutate(ageatV3 = (v3_date - dob),
  #        ageatlastvac = last_primary_vac - dob) %>%
  # select(-dob) %>%
  # add same end date for now (should be updated from Paul)
  mutate(fu_end_date = ymd('2020-03-31'))  


# Get total children per arm
total_children <- children_raw %>%
  group_by(arm) %>%
  summarise(n_children = n_distinct(rid), .groups = "drop")

primary <- primary_raw %>%
  left_join(children) %>%
  mutate(age_at_inf = dcontact - dob) %>%
  group_by(arm) %>%
  mutate(mean_age_at_inf = mean(age_at_inf, na.rm = TRUE),
         week = tsibble::yearweek(dcontact)) %>% ungroup() %>%
  # Get first infection per child 
  arrange(rid, dcontact) %>%
  group_by(rid) %>%
  mutate(first_inf = ifelse(poutcome == 1 & row_number() == min(which(poutcome == 1)),1, 0)) %>%
  ungroup() %>%
  left_join(total_children) %>%
  group_by(arm) %>%
  arrange(dcontact, arm) %>%
  mutate(cum_cases = cumsum(poutcome),
         dailycum_inci = cum_cases / n_children) %>%
  # Calculate monthly incidence
  mutate(month = lubridate::month(dcontact),
         year = lubridate::year(dcontact),
         monthyear = paste0(month, "-", year)) %>%
  group_by(month, year, arm) %>%
  mutate(cases_per_month = sum(poutcome))


delivery <- delivery_raw %>% # has dob
  left_join(children)

# has age_months which is age in months at v1_date
delivery_detail <- delivery_detail_raw %>%
  left_join(children) 

serology <- serology_raw %>%
  left_join(children) %>%
  mutate(prepost_v3 = ifelse(sdate <= v3_date, 'pre', 
                             ifelse(sdate > v3_date, 'post', NA)),
         time_since_vac = sdate - last_primary_vac)

weekly_raw <- rbind(weekly1_raw, weekly2_raw, weekly3_raw) 
weekly <- weekly_raw %>%
  left_join(children) %>%
  mutate(time_since_vac = ifelse(pf_asex_fresult == 1 & arm != 'smc', dateweekly - v3_date, NA),
         time_since_smc = ifelse(pf_asex_fresult == 1 & arm != 'rtss', dateweekly - y1p1d1_date_received, NA))


# Analysis dataset
# information about every child in study, primary outcomes, and delivery details 
# mitt - which is anyone who got th first vaccine
mitt <- children %>%
  left_join(delivery_detail_raw, by = c('arm','sex','country','rid','v1_date','v2_date','v3_date','last_primary_vac')) %>%
  left_join(primary_raw, by = c('rid')) %>%
  # Do same cleaning as above for children
  # Get approx date of birth
  mutate(dob = v1_date - age_months_v1 * 30)%>%
  # Get age at 3rd primary series vaccination and if not the 3rd, then the most recent one
  mutate(ageatV3 = (v3_date - dob),
         ageatlastvac = last_primary_vac - dob) %>%
  # and for primary
  # get age at each infection in days 
  mutate(age_at_inf = dcontact - dob) %>%
  # find average age at infection by arm and create a week variable 
  group_by(arm) %>%
  mutate(mean_age_at_inf = mean(age_at_inf, na.rm = TRUE),
         week = tsibble::yearweek(dcontact)) %>% ungroup() %>%
  # Get first infection per child 
  arrange(rid, dcontact) %>%
  group_by(rid) %>%
  mutate(first_inf = if (any(poutcome == 1, na.rm = TRUE)) {
    as.integer(poutcome == 1 & row_number() == min(which(poutcome == 1)))
  } else {
    0L
  }) %>%
  ungroup() %>%
  # join with total children df to get #kids in each arm 
  left_join(total_children) %>%
  group_by(arm) %>%
  # calculate cumulative cases per arm and also incidence per person not per PY
  arrange(dcontact, arm) %>%
  mutate(cum_cases = cumsum(poutcome),
         dailycum_inci = cum_cases / n_children) %>%
  # Calculate monthly incidence
  mutate(month = lubridate::month(dcontact),
         year = lubridate::year(dcontact),
         monthyear = paste0(month, "-", year)) %>%
  group_by(month, year, arm) %>%
  mutate(cases_per_month = sum(poutcome)) %>%
  ungroup() %>%
  group_by(rid, month, year, monthyear) %>%
  mutate(n_infections_per_month_pp = n())


# Cleaning the person time file 
primary_persontime <- primary_persontime %>%
  mutate(person_years = X_t - X_t0) %>%
  mutate(arm = case_when(
    arm == 1 ~ 'smc',
    arm == 2 ~ 'rtss',
    arm == 3 ~ 'both',
    TRUE ~ NA),
    country = case_when(
      country == 1 ~ 'BF',
      country == 2 ~ 'Mali',
      TRUE ~ NA),
    arm_smcref = factor(arm, levels = c('smc','rtss','both')),
    arm_rtssref = factor(arm, levels = c('rtss','smc','both'))) %>%
  rename(end_time = X_t,
         start_time = X_t0,
         event = X_d) 


dir.create('data/')
cyphr::encrypt(saveRDS(primary_persontime, file = 'data/primary_persontime.rds'), key)
cyphr::encrypt(saveRDS(mitt, file = 'data/mitt.rds'), key)
cyphr::encrypt(saveRDS(weekly, file = 'data/weekly.rds'), key)
cyphr::encrypt(saveRDS(serology, file = 'data/serology.rds'), key)
cyphr::encrypt(saveRDS(delivery_detail, file = 'data/delivery_detail.rds'), key)
cyphr::encrypt(saveRDS(primary, file = 'data/primary.rds'), key)
cyphr::encrypt(saveRDS(children, file = 'data/children.rds'), key)
