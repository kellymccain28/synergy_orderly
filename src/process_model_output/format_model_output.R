format_model_output <- function(model_data, 
                                start_fu = as.Date('2017-04-01')){
  
  # First, convert to date format adnd get vaccination date 
  model_data <- model_data %>%
    
    # make date vars be Date format 
    mutate(across(c(time_ext, infectious_bite_day, BSinfection_day, detection_day),
                  ~ as.Date(.x, origin = '2017-04-01'))) %>%
    
    # add date for vaccination 
    mutate(vaccinate_date = start_fu + vaccination_day) %>%# vaccination day is negative if before the start of follow-up
    
    # Find year of infection 
    mutate(infection_year = case_when(
      detection_day >= start_fu & detection_day < start_fu + 365 ~ 1,
      detection_day >= start_fu + 365 & detection_day < start_fu + 365*2 ~ 2, 
      detection_day >= start_fu + 365*2 & detection_day < start_fu + 365*3 ~ 3,
      TRUE ~ NA
    ))
  
  # Function to format rows with detectable cases specific years-- 
  format_cases <- function(year){
    
    censor_date = as.Date(start_fu + 365*year)
    
    df <- model_data %>%
      
      # First, filter to only include rows with infections 
      filter(#!is.na(detectable) & detectable == 1 &
        ( detection_day <= censor_date & 
            detection_day > (censor_date - 365) ) #| 
        # (BSinfection_day <= censor_date &
        #    BSinfection_day > (censor_date - 365) &
        #    is.na(detection_day)) 
      ) %>%
      
      # Group by child and arrange by infection time
      group_by(child_id) %>%
      arrange(child_id, detection_day) %>%
      
      # Create recurrent event structure
      mutate(
        
        # Calculate start time for each interval
        start_date = case_when(
          row_number() == 1 & year == 1 ~ start_fu, # First infection starts at origin
          row_number() == 1 & year != 1 ~ as.Date(start_fu + 365 * (year - 1)),
          TRUE ~ lag(detection_day)#time_length(interval(start_fu, lag(detection_day)), unit = 'years') #(lag(detection_day) - start_fu) / 365.25
        ),
        
        # End time is when this infection was detected or the end of the year, whichever is first 
        end_date = if_else(detection_day < censor_date, detection_day, censor_date),#time_length(interval(start_fu, detection_day), unit = 'years') ,#(detection_day - start_fu) / 365.25,
        
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
      group_by(child_id) %>%
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
        censoring_row$infectious_bite_day <- NA
        censoring_row$time_ext <- NA
        censoring_row$detection_day <- NA
        censoring_row$t_toreach_threshold <- NA
        censoring_row$detectable <- 0
        censoring_row$infection_year <- NA
        
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
  
  # Get all unique child ids 
  all_children <- unique(metadata_child$child_id)
  
  # Year framework 
  year_framework <- tibble(
    child_id = rep(all_children, each = 3),
    year = rep(1:3, times = length(all_children))
  ) %>%
    mutate(year_start = start_fu + 365 * (year - 1),
           year_end = start_fu + 365 * year) %>%
    left_join(metadata_child)
  
  # Join with exsiting combinations 
  existing_combinations <- all_cases %>%
    distinct(child_id, year)
  
  # Find the missing ones 
  missing_combinations <- year_framework %>%
    anti_join(existing_combinations, by = c('child_id','year'))
  
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
    left_join(model_data %>% distinct(intervention, children_in_group))
  
  # Bind all datasets together 
  person_time <- rbind(all_cases, censoring_intervals) %>%
    arrange(child_id) %>%
    # Get start and end times in years since origin
    mutate(origin = start_fu, 
           start_time = time_length(interval(origin, start_date), unit = 'years'),
           end_time = time_length(interval(origin, end_date), unit = 'years'),
           person_years = end_time - start_time)
  
  return(person_time)
}