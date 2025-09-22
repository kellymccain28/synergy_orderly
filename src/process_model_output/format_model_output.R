format_model_output <- function(model_data, 
                                cohort, # 'generic' or 'trial' -- this influences the start of follow-up
                                start_cohort = as.Date('2017-04-01'),
                                simulation){
  
  # First, convert to date format adnd get vaccination date 
  model_data <- model_data %>%
    
    filter(sim_id == simulation) %>%
    
    # make date vars be Date format 
    mutate(across(c(time_ext, infectious_bite_day, BSinfection_day, detection_day),
                  ~ as.Date(.x, origin = start_cohort))) %>%
    
    # add date for vaccination 
    mutate(vaccinate_date = start_cohort + vaccination_day) %>%# vaccination day is negative if before the start of follow-up
    
    # Find year of infection 
    mutate(infection_year = case_when(
      detection_day >= start_cohort & detection_day < start_cohort + 365 ~ 1,
      detection_day >= start_cohort + 365 & detection_day < start_cohort + 365*2 ~ 2, 
      detection_day >= start_cohort + 365*2 & detection_day < start_cohort + 365*3 ~ 3,
      TRUE ~ NA
    )) %>%
    #Filter to remove infections that were detected prior to vaccination (which is first day of follow-up)
    filter(v1_date < detection_day | is.na(detection_day))
    
    
  
  # Function to format rows with detectable cases specific years-- 
  format_cases <- function(year){
    
    censor_date = as.Date(start_cohort + 365*year)
    
    df <- model_data %>%
      
      # First, filter to only include rows with infections 
      filter(
         detection_day <= censor_date & 
            detection_day > (censor_date - 365) 
      ) %>%
      
      # Group by child and arrange by infection time
      group_by(rid) %>%
      arrange(rid, detection_day) %>%
      
      # Create recurrent event structure
      mutate(

        # if the cohort argument is 'generic' then use the satart of the cohort for everyone, otherwise, it is the first vaccination day
        start_fu_date = case_when(
          cohort == 'generic' ~ start_cohort,
          TRUE ~ v1_date),
        
        # Calculate start time for each interval
        start_date = case_when(
          row_number() == 1 & year == 1 ~ as.Date(start_fu_date), # First infection starts at origin
          row_number() == 1 & year != 1 ~ as.Date(start_cohort + 365 * (year - 1)), # start dates for following years are on these same days for everyone
          TRUE ~ lag(detection_day)#time_length(interval(start_cohort, lag(detection_day)), unit = 'years') #(lag(detection_day) - start_cohort) / 365.25
        ),
        
        # End time is when this infection was detected or the end of the year, whichever is first 
        end_date = if_else(detection_day < censor_date, detection_day, censor_date),#time_length(interval(start_cohort, detection_day), unit = 'years') ,#(detection_day - start_cohort) / 365.25,
        
        # All these are events (infections)
        event = 1,
        
        # year number
        year = year,
        
        # Use detection_day as contact date
        dcontact = if_else(detection_day < censor_date, detection_day, censor_date),
        
        poutcome = 1 # All are infections
      ) %>%
      
      ungroup() %>%
      
      # Add censoring intervals for people who had infections
      # (interval from last infection to end of follow-up) -- so during this interval there is no infection so all infection-related values should be 0 
      group_by(rid) %>%
      group_modify(~{
        last_infection <- .x[nrow(.x), ]
        
        # Create censoring row
        censoring_row <- last_infection
        censoring_row$start_date <- as.Date(last_infection$end_date)
        censoring_row$end_date <- censor_date
        censoring_row$event <- 0
        # censoring_row$year <- nrow(.x) + 1
        censoring_row$dcontact <- censoring_row$start_date 
        censoring_row$poutcome <- NA_real_
        censoring_row$threshold_day <- NA
        censoring_row$BSinfection_day <- NA
        censoring_row$infectious_bite_day <- as.Date(NA)
        censoring_row$time_ext <- as.Date(NA)
        censoring_row$detection_day <- NA
        censoring_row$t_toreach_threshold <- NA
        censoring_row$detectable <- 0
        censoring_row$infection_year <- NA
        censoring_row$start_fu_date <- censoring_row$v1_date
        
        # Combine infections + censoring
        bind_rows(.x, censoring_row)
      })
  }
  
  # Format cases by year 
  cases_yr1 <- format_cases(year = 1)
  cases_yr2 <- format_cases(year = 2)
  cases_yr3 <- format_cases(year = 3)
  
  # Bind all cases together
  all_cases <- bind_rows(cases_yr1, cases_yr2, cases_yr3)
  
  # Filter metadata_child to be only for simulation sim
  metadata_child <- metadata_child %>%
    filter(sim_id == simulation)
  
  # Get all unique child ids 
  all_children <- unique(metadata_child$rid)
  
  # Year framework 
  year_framework <- tibble(
    rid = rep(all_children, each = 3),
    year = rep(1:3, times = length(all_children))
  ) %>%
    mutate(year_start = start_cohort + 365 * (year - 1),
           year_end = start_cohort + 365 * year) %>%
    left_join(metadata_child)
  
  # Join with exsiting combinations 
  existing_combinations <- all_cases %>%
    distinct(rid, year)
  
  # Find the missing ones 
  missing_combinations <- year_framework %>%
    anti_join(existing_combinations, by = c('rid','year'))
  
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
      infection_year = NA
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
  
  message('finished ', simulation)
  
  return(person_time)
}