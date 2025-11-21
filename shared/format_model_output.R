# Function to format the model output (infection_records) to have person time 
format_model_output <- function(model_data, 
                                cohort, # 'generic' or 'trial' -- this influences the start of follow-up
                                start_cohort = as.Date('2017-04-01'),
                                simulation){
  
  # First, convert to date format and get vaccination date 
  model_data <- model_data %>%
    
    filter(sim_id == simulation) %>%
    
    # make date vars be Date format 
    mutate(across(c(time_ext, infectious_bite_day, BSinfection_day, detection_day),
                  ~ as.Date(.x, origin = start_cohort))) %>%
    
    # before anything else, remove any bites that didn't turn into infections
    # and those detected before vaccination (first day of follow-up)
    filter(!is.na(detection_day) & v1_date <= detection_day) %>%
    
    # add date for vaccination 
    mutate(vaccinate_date = start_cohort + vaccination_day) %>% # vaccination day is negative if before the start of follow-up
    
    # Find year of infection 
    mutate(infection_year = floor(as.numeric(detection_day - start_cohort) / 365) + 1
  )
  
  # Function to format rows with detectable cases by specific years-- 
  format_cases_by_year <- function(model_data, yr){
    
    censor_date = as.Date(start_cohort + 365*yr)
    
    df <- model_data %>%
      
      # First, filter to only include rows with infections 
      filter(
         detection_day <= censor_date & 
            detection_day > (censor_date - 365) 
      ) %>%
      
      # Group by child and arrange by infection time
      group_by(rid) %>%
      arrange(detection_day) %>%
      
      # Create recurrent event structure
      mutate(

        # if the cohort argument is 'generic' then use the start of the cohort for everyone, otherwise, it is the first vaccination day
        start_fu_date = case_when(
          cohort == 'generic' ~ start_cohort,
          TRUE ~ v1_date),
        
        # Calculate start time for each interval
        start_date = case_when(
          row_number() == 1 & yr == 1 ~ as.Date(start_fu_date), # First infection starts at origin
          row_number() == 1 & yr != 1 ~ as.Date(start_cohort + 365 * (yr - 1)), # start dates for following years are on these same days for everyone
          TRUE ~ lag(detection_day)#time_length(interval(start_cohort, lag(detection_day)), unit = 'years') #(lag(detection_day) - start_cohort) / 365.25
        ),
        
        # End time is when this infection was detected or the end of the year, whichever is first 
        end_date = pmin(detection_day, censor_date),#if_else(detection_day < censor_date, detection_day, censor_date),#time_length(interval(start_cohort, detection_day), unit = 'years') ,#(detection_day - start_cohort) / 365.25,
        
        # All these are events (infections)
        event = 1,
        
        # year number
        year = yr,
        
        # Use detection_day as contact date
        dcontact = end_date,#if_else(detection_day < censor_date, detection_day, censor_date),
        
        poutcome = 1 # All are infections
      ) %>%
      
      ungroup() #%>%
    
    last_infections <- df %>%
      group_by(rid) %>%
      slice_max(end_date, n = 1, with_ties = FALSE) %>% # shoudl this be end or start??
      ungroup()

    # Create censoring rows
    censoring_rows <- last_infections %>%
      mutate(
        start_date = end_date,
        end_date = censor_date,
        event = 0,
        dcontact = start_date,
        poutcome = NA_real_,
        threshold_day = NA,
        BSinfection_day = as.Date(NA),
        infectious_bite_day = as.Date(NA),
        time_ext = as.Date(NA),
        detection_day = as.Date(NA),
        t_toreach_threshold = NA,
        detectable = 0,
        infection_year = NA,
        start_fu_date = v1_date
      )

    # Combine infections + censoring
    bind_rows(df, censoring_rows) %>%
      arrange(rid, start_date)

  }
  
  # Format cases by year 
  # cases_yr1 <- format_cases_by_year(year = 1)
  # cases_yr2 <- format_cases(year = 2)
  # cases_yr3 <- format_cases(year = 3)
  
  # Bind all cases together
  all_cases <- map_dfr(1:3, ~format_cases_by_year(model_data, .x))#bind_rows(cases_yr1, cases_yr2, cases_yr3)
  
  # Get all unique child ids 
  all_children <- unique(metadata_df$rid)
  
  # Children without infections
  # Year framework 
  year_framework <- tibble(
    rid = rep(all_children, each = 3),
    year = rep(1:3, times = length(all_children))
  ) %>%
    mutate(year_start = start_cohort + 365 * (year - 1),
           year_end = start_cohort + 365 * year) %>%
    left_join(metadata_df) %>%
    left_join(all_cases %>% distinct(rid, smckillvec, smc_dose_days, v1_date))
  
  # Join with exsiting combinations and find the missing ones 
  missing_combinations <- year_framework %>%
    anti_join(all_cases %>%
                distinct(rid, year), by = c('rid','year'))
  
  # Create censoring intervals for missing combinations 
  censoring_intervals <- missing_combinations %>%
    mutate(
      start_date = year_start, 
      end_date = year_end, 
      event = 0, 
      dcontact = year_start, 
      poutcome = NA, 
      time_ext = NA, 
      infectious_bite_day = NA, 
      BSinfection_day = NA, 
      threshold_day = NA, 
      detection_day = NA, 
      t_toreach_threshold = NA, 
      detectable = 0, 
      vaccinate_date = as.Date(vaccination_day, origin = as.Date('2017-04-01')),
      infection_year = NA,
      prob_bite = NA, 
      recovery_day = NA, 
      cumul_inf = NA, 
      start_fu_date = case_when(
        cohort == 'generic' ~ start_cohort,
        TRUE ~ v1_date),
      t = NA,
      sim_id = simulation,
      treatment_day = NA,
      receives_treatment = NA,
      treatment_efficacy = NA, 
      treatment_successful = NA
    ) %>% 
    select(-year_start, -year_end) %>%
    # add children_in_group var
    left_join(model_data %>% distinct(arm, children_in_group))
  
  # Bind all datasets together 
  person_time <- rbind(all_cases, censoring_intervals) %>%
    arrange(rid) %>%
    # Get start and end times in years since origin
    mutate(origin = start_cohort, 
           start_time = time_length(interval(origin, start_date), unit = 'years'),
           end_time = time_length(interval(origin, end_date), unit = 'years'),
           person_years = end_time - start_time)
  
  # saveRDS(person_time, file = paste0('outputs/processed_', simulation, ".rds"))
  message('finished ', simulation)
  
  return(person_time)
}