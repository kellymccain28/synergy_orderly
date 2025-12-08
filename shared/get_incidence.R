# Functions to get incidnece per 1000 person months 
# output will have cases per month, person months per month, and cases/1000 person months
# can be run on model output or on trial data

get_incidence <- function(model = TRUE, 
                          df_children, #if model, then `metadata_child`, if trial then `children`
                          casedata #if model, then `model_output`, if trial then `mitt`
                          
){
  # Define function to expand to get all the different months in the dataset
  # make_child_months <- function(rid, arm, start_date, fu_end_date) {
  #   month_seq <- seq(floor_date(start_date, "month"), floor_date(fu_end_date, "month"), by = "month")
  #   tibble(
  #     rid = rid,
  #     arm = arm,
  #     month = floor_date(month_seq, "month")
  #   )
  # }
  
  # expand all children 
  person_months_df <- df_children %>%
    mutate(fu_end_date = ymd('2020-03-31'),# this should change to when the follow up actually ended
           start_date = v1_date) %>%
    
    ungroup() %>% 
    mutate(month = map2(start_date, fu_end_date,
                        ~seq(floor_date(.x, "month"),
                             floor_date(.y, "month"),
                             by = 'month'))) %>%
    unnest(month)

  # person_months_df <- person_months_df %>%
  #   ungroup() %>%
  #   select(rid, arm, start_date, fu_end_date) %>%
  #   pmap(~ make_child_months(..1, ..2, as.Date(..3), as.Date(..4))) %>%
  #   bind_rows()
  # 
  
  # person_months_df2 <- person_months_df %>%
  #   left_join(df_children %>% 
  #               select(rid, country), by = c('country','rid'))
  
  # Calculate actual person-time per child-month
  # for each month, we need to find max an dmin -- so if the start date is=or > the floor date of the month and < last day of that month, 
  # then we take either the first day of the month or the start date. if the start date is > last da of taht month, then we take the 
  # last day of that month (I think????)
  # for the end date, i think it's fine 
  person_months_df3 <- person_months_df %>%
    mutate(
      month_start = pmax(floor_date(month, "month"), v1_date),
      month_end = pmin(ceiling_date(month, "month") - days(1), fu_end_date),
      days_at_risk = time_length(interval(month_start, month_end), unit = 'days'),
      person_months = time_length(interval(month_start, month_end), unit = 'months'),
      person_years = time_length(interval(month_start, month_end), unit = 'years'),
    )
  
  # Summarize by calendar month and arm
  persontime_bymonth <- person_months_df3 %>%
    mutate(year = lubridate::year(month), 
           month_num = month(month),
           yearmonth = zoo::as.yearmon(month)) %>%
    group_by(year, month, yearmonth, arm) %>%
    summarise(
      person_months = sum(person_months),
      person_years = sum(person_years),
      .groups = "drop"
    )
  
  if(model){
    # Calculate incidence per month ----
    cases_by_month <- casedata %>%
      filter(poutcome == 1) %>%
      mutate(month = floor_date(detection_day, "month"), 
             month_num = lubridate::month(detection_day),
             year = lubridate::year(detection_day),
             yearmonth = zoo::as.yearmon(month)) %>%
      group_by(arm, month, month_num, year, yearmonth) %>%
      summarize(n_cases = sum(poutcome), .groups = 'drop') 
  } else {
    # Calculate incidence per month ----
    cases_by_month <- casedata %>%
      filter(poutcome == 1 )%>%
      mutate(month_num = month, 
             yearmonth = zoo::as.yearmon(dcontact),
             month = floor_date(dcontact, "month"),
             year = lubridate::year(dcontact)) %>%
      group_by(arm, month_num, month, year, yearmonth) %>%
      summarize(n_cases = sum(poutcome), .groups = 'drop')
  }
  
  # join t person time data 
  monthly_inci <- persontime_bymonth %>%
    left_join(cases_by_month, by = c('year','month','yearmonth','arm')) %>%
    mutate(
      n_cases = replace_na(n_cases, 0),
      incidence_per_1000pm = (n_cases / person_months) * 1000,
      rate = n_cases / person_months, 
      se = sqrt(n_cases) / person_months,
      lower = (qchisq(0.025, 2 * n_cases) / 2) / person_months, 
      upper = (qchisq(0.975, 2 * n_cases + 2) / 2) / person_months,
      incidence_per_1000pm = rate * 1000,
      lower_per_1000 = lower * 1000,
      upper_per_1000 = upper * 1000,
    ) %>%
    mutate(date = make_date(year, month_num, 1))
  
  return(monthly_inci)
  
}


get_incidence_simpler <- function(model_data, 
                                  metadata_df, 
                                  start_cohort = as.Date('2017-04-01'), 
                                  time_horizon) {#weekly or monthly
  
  child_counts <- metadata_df %>%
    count(arm, country, name = 'pop')#, sim_id
  
  # Get case counts by month 
  case_counts <- model_data %>%
    # mutate(detection_day = as.Date(detection_day)) %>% # for trial data
    mutate(detection_day = as.Date(detection_day, origin= start_cohort)) %>%
    filter(!is.na(detection_day)) %>%
    filter(detection_day >= start_cohort)  %>%
    filter(detection_day <= (start_cohort + 365*3)) %>% 
    mutate(yearmonth = floor_date(detection_day, "month"),
           week = floor_date(detection_day, "week")) 
  
  if(time_horizon == 'monthly'){
    output <- case_counts %>% 
      count(arm, country, yearmonth,  name = "n_cases") %>% #sim_id,
      tidyr::complete(yearmonth = seq(min(yearmonth), max(yearmonth), by = '1 month'), 
                      fill = list(cases = 0,
                                  inci = 0))  %>%
      left_join(child_counts) %>%
      mutate( 
        # Approximate person-months (cohort size × 1 month) 
        person_months = pop,
        incidence_per_1000pm = (n_cases / person_months) * 100 )
    
  }
  
  if(time_horizon == 'weekly'){
    output <- case_counts %>%
      count(arm, country, week,  name = "n_cases") %>% #sim_id,
      tidyr::complete(week = seq(min(week), max(week), by = '1 week'), 
                      fill = list(cases = 0,
                                  inci = 0)) %>%
      left_join(child_counts) %>%
      mutate( 
        # Approximate person-months (cohort size × 1 month) 
        person_weeks = pop,
        incidence_per_1000pw = (n_cases / person_weeks) * 100 )
  }
  return(output) 
}


# get_poisson_inci <- function(df_children = children,
#                              casedata = mitt){
#   # can use a poisson regression to get monthly incidence but need to format the data differently than how it is formatted now 
#   
#   # expand all children 
#   person_months_df <- df_children %>%
#     mutate(fu_end_date = ymd('2020-03-31'),# this should change to when the follow up actually ended
#            start_date = v1_date) %>%
#     
#     ungroup() %>% 
#     mutate(month = map2(start_date, fu_end_date,
#                         ~seq(floor_date(.x, "month"),
#                              floor_date(.y, "month"),
#                              by = 'month'))) %>%
#     unnest(month)
#   
#   # person_months_df2 <- person_months_df %>%
#   #   left_join(df_children %>% 
#   #               select(rid, country), by = 'rid')
#   
#   # Calculate actual person-time per child-month
#   # for each month, we need to find max an dmin -- so if the start date is=or > the floor date of the month and < last day of that month, 
#   # then we take either the first day of the month or the start date. if the start date is > last da of taht month, then we take the 
#   # last day of that month (I think????)
#   # for the end date, i think it's fine 
#   person_months_df3 <- person_months_df %>%
#     mutate(
#       month_start = pmax(floor_date(month, "month"), v1_date),
#       month_end = pmin(ceiling_date(month, "month") - days(1), fu_end_date),
#       days_at_risk = time_length(interval(month_start, month_end), unit = 'days'),
#       person_months = time_length(interval(month_start, month_end), unit = 'months'),
#       person_years = time_length(interval(month_start, month_end), unit = 'years'),
#     )
#   
#   # Summarize by calendar month and arm
#   persontime_bymonth <- person_months_df3 %>%
#     mutate(year = lubridate::year(month), 
#            month_num = month(month),
#            yearmonth = zoo::as.yearmon(month)) %>%
#     group_by(year, month, yearmonth, arm) %>%
#     summarise(
#       person_months = sum(person_months),
#       person_years = sum(person_years),
#       .groups = "drop"
#     )
#   
#   # Calculate incidence per month ----
#   cases_by_month <- casedata %>%
#     filter(poutcome == 1 )%>%
#     mutate(month_num = month, 
#            yearmonth = zoo::as.yearmon(dcontact),
#            month = floor_date(dcontact, "month"),
#            year = lubridate::year(dcontact)) %>%
#     group_by(arm, month_num, month, year, yearmonth) %>%
#     summarize(n_cases = sum(poutcome), .groups = 'drop')
#   
#   monthly_inci <- persontime_bymonth %>%
#     left_join(cases_by_month, by = c('year','month','yearmonth','arm')) 
#   
#   fit <- glm(
#     n_cases ~ arm + yearmonth,                  # predictors
#     offset = log(person_months),            # exposure (person-time)
#     family = poisson,                       # Poisson model for count data
#     data = monthly_inci
#   )
#   
#   # summary(fit)
#   newdat <- monthly_inci %>%
#     distinct(arm, yearmonth, person_months) %>%
#     group_by(arm, yearmonth) %>%
#     summarise(person_months = mean(person_months), .groups = "drop")
#   pred <- predict(fit, newdata = newdat, type = "link", se.fit = TRUE)
#   
#   newdat <- newdat %>%
#     mutate(
#       fit = pred$fit,
#       se = pred$se.fit,
#       lower = fit - 1.96 * se,
#       upper = fit + 1.96 * se,
#       # back-transform from log scale to rate per person-month
#       rate = exp(fit),
#       lower_rate = exp(lower),
#       upper_rate = exp(upper),
#       # convert to incidence per 1,000 person-months
#       incidence_per_1000 = rate * 1000,
#       lower_per_1000 = lower_rate * 1000,
#       upper_per_1000 = upper_rate * 1000
#     )
#   
#   ggplot(newdat, aes(x = yearmonth, y = rate, color = arm)) +
#     geom_line(size = 1) +
#     geom_ribbon(aes(ymin = lower_rate, ymax = upper_rate, fill = arm),
#                 alpha = 0.2, color = NA) +
#     labs(
#       x = "Month",
#       y = "Monthly malaria incidence rate\n(per 1000 person-months)",
#       title = "Estimated monthly incidence by trial arm",
#       color = "Trial arm", fill = "Trial arm"
#     ) +
#     theme_minimal(base_size = 13)
#   
# }