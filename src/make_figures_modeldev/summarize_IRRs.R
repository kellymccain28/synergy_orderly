# Function to use in 'plot_1-IRR_average.R' which uses time aggregation of 'year' or 'halfyear', or in 
# 'plot_1-IRR.R' which uses 'yearmonth' 
# This function will aggregate incidence as needed, calculate IRRs, then summarize by the time unit 
# to get median and 95% CrIs of the IRRs, efficacy (1-IRR), and the ratio of efficacy between expected and model-predicted using bootstrapping

summarize_IRRs <- function(outputsfolder, 
                           agg_unit){
  library(purrr)
  path <- paste0('R:/Kelly/synergy_orderly/src/', cohort_folder, '/outputs/')
  inci <- readRDS(paste0(path, outputsfolder, '/incidence.rds'))
  
  if(agg_unit == 'year'){
    # Get annual and overall incidence 
    inci_annual <- inci %>%
      filter(!is.na(date)) %>%
      mutate(time_value = case_when(date < '2018-04-01' & date > '2017-05-01' ~ 'Jun 2017-Mar 2018',
                                   date < '2019-04-01' ~ 'Apr 2018-Mar 2019',
                                   date < '2020-04-01' ~ 'Apr 2019-Mar 2020'),
             time_value_num = case_when(date < '2018-04-01' & date > '2017-05-01' ~ '1',
                                        date < '2019-04-01' ~ '2',
                                        date < '2020-04-01' ~ '3')) %>%
      group_by(time_value, time_value_num, arm, sim_id) %>%
      summarize(person_months = sum(person_months),
                n_cases = sum(n_cases)) %>%
      mutate(incidence_per_1000pm = n_cases / person_months * 1000,
             time_value = as.character(time_value),
             time_unit = agg_unit)
    
    inci_overall <- inci %>%
      group_by(arm, sim_id) %>%
      summarize(person_months = sum(person_months),
                n_cases = sum(n_cases)) %>%
      mutate(incidence_per_1000pm = n_cases / person_months * 1000,
             time_value = 'Overall',
             time_value_num = 'Overall',
             time_unit = agg_unit)
    
    inci <- rbind(inci_annual, inci_overall) %>%
      mutate(time_value = factor(time_value, levels = c('Jun 2017-Mar 2018',
                                                        'Apr 2018-Mar 2019',
                                                        'Apr 2019-Mar 2020',
                                                        'Overall')))
    
  } else if (agg_unit == 'halfyear'){
    inci_annual <- inci %>%
      filter(!is.na(date)) %>%
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
      group_by(time_value, time_value_num, arm, sim_id) %>%
      summarize(person_months = sum(person_months),
                n_cases = sum(n_cases)) %>%
      mutate(incidence_per_1000pm = n_cases / person_months * 1000,
             time_value = as.character(time_value),
             time_unit = agg_unit)
    
    inci_overall <- inci %>%
      group_by(arm, sim_id) %>%
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
  } else if( agg_unit == 'yearmonth'){
    inci <- inci %>%
      mutate(time_value = yearmonth, 
             time_value_num = yearmonth,
             time_unit = agg_unit)
  }
  
  # Pivot wider 
  inci_wide <- inci %>%
    split(.$sim_id) %>%
    map_dfr(~ .x %>%
              select(arm, time_value, time_value_num, time_unit, 
                     person_months, incidence_per_1000pm) %>%
              pivot_wider(
                names_from = arm,
                values_from = c(person_months, incidence_per_1000pm),
                id_cols = c(time_value, time_value_num, time_unit)
              ),
            .id = "sim_id")
  
  # Get IRRs
  inci_wide <- inci_wide %>%
    mutate(rtss_none_irr = (incidence_per_1000pm_rtss / incidence_per_1000pm_none),
           smc_none_irr = (incidence_per_1000pm_smc / incidence_per_1000pm_none),
           both_none_irr = (incidence_per_1000pm_both / incidence_per_1000pm_none),
           rtss_smc_irr = (incidence_per_1000pm_rtss / incidence_per_1000pm_smc),
           both_smc_irr = (incidence_per_1000pm_both / incidence_per_1000pm_smc),
           both_rtss_irr = (incidence_per_1000pm_both / incidence_per_1000pm_rtss),
           smc_rtss_irr = (incidence_per_1000pm_smc / incidence_per_1000pm_rtss)  )%>%
    mutate(expected_efficacy = 1 - (rtss_none_irr * smc_none_irr),
           ratio_pred_exp = (1-both_none_irr) / expected_efficacy)
  
  # Calculate median and IQR for each metric by time aggregation unit 
  inci_summary <- inci_wide %>%
    group_by(time_value, time_value_num, time_unit) %>%
    reframe(
      # For incidence, can take direct median and quantiles
      # incidence smc 
      incidence_smc_median = median(incidence_per_1000pm_smc, na.rm = TRUE),
      incidence_smc_q25 = quantile(incidence_per_1000pm_smc, 0.25, na.rm = TRUE),
      incidence_smc_q75 = quantile(incidence_per_1000pm_smc, 0.75, na.rm = TRUE),
      incidence_smc_q025 = quantile(incidence_per_1000pm_smc, 0.025, na.rm = TRUE),
      incidence_smc_q975 = quantile(incidence_per_1000pm_smc, 0.975, na.rm = TRUE),
      # incidence rtss 
      incidence_rtss_median = median(incidence_per_1000pm_rtss, na.rm = TRUE),
      incidence_rtss_q25 = quantile(incidence_per_1000pm_rtss, 0.25, na.rm = TRUE),
      incidence_rtss_q75 = quantile(incidence_per_1000pm_rtss, 0.75, na.rm = TRUE),
      incidence_rtss_q025 = quantile(incidence_per_1000pm_rtss, 0.025, na.rm = TRUE),
      incidence_rtss_q975 = quantile(incidence_per_1000pm_rtss, 0.975, na.rm = TRUE),
      # incidence none 
      incidence_none_median = median(incidence_per_1000pm_none, na.rm = TRUE),
      incidence_none_q25 = quantile(incidence_per_1000pm_none, 0.25, na.rm = TRUE),
      incidence_none_q75 = quantile(incidence_per_1000pm_none, 0.75, na.rm = TRUE),
      incidence_none_q025 = quantile(incidence_per_1000pm_none, 0.025, na.rm = TRUE),
      incidence_none_q975 = quantile(incidence_per_1000pm_none, 0.975, na.rm = TRUE),
      # incidence both 
      incidence_both_median = median(incidence_per_1000pm_both, na.rm = TRUE),
      incidence_both_q25 = quantile(incidence_per_1000pm_both, 0.25, na.rm = TRUE),
      incidence_both_q75 = quantile(incidence_per_1000pm_both, 0.75, na.rm = TRUE),
      incidence_both_q025 = quantile(incidence_per_1000pm_both, 0.025, na.rm = TRUE),
      incidence_both_q975 = quantile(incidence_per_1000pm_both, 0.975, na.rm = TRUE),
      
      # Bootstrap for IRRs converted to efficacy
      both_smc = list(bootstrap_metric(both_smc_irr, transform_fn = function(x) 1 - x)),
      rtss_none = list(bootstrap_metric(rtss_none_irr, transform_fn = function(x) 1 - x)),
      both_rtss = list(bootstrap_metric(both_rtss_irr, transform_fn = function(x) 1 - x)),
      smc_none = list(bootstrap_metric(smc_none_irr, transform_fn = function(x) 1 - x)),
      both_none = list(bootstrap_metric(both_none_irr, transform_fn = function(x) 1 - x)),
      rtss_smc = list(bootstrap_metric(rtss_smc_irr, transform_fn = function(x) 1 - x)),
      smc_rtss = list(bootstrap_metric(smc_rtss_irr, transform_fn = function(x) 1 - x)),
      # Bootstrap for expected efficacy
      expected_efficacy = list(bootstrap_metric(expected_efficacy)),
      
      # Bootstrap for ratio
      ratio_pred_exp = list(bootstrap_metric(ratio_pred_exp))
    ) %>%
    # Unpack all the bootstrap results
    mutate(
      # Both vs SMC
      both_smc_median = sapply(both_smc, `[`, 1),
      both_smc_q025 = sapply(both_smc, `[`, 2),
      both_smc_q975 = sapply(both_smc, `[`, 3),
      both_smc_q25 = sapply(both_smc, `[`, 4),
      both_smc_q75 = sapply(both_smc, `[`, 5),
      
      # RTSS vs none
      rtss_none_median = sapply(rtss_none, `[`, 1),
      rtss_none_q025 = sapply(rtss_none, `[`, 2),
      rtss_none_q975 = sapply(rtss_none, `[`, 3),
      rtss_none_q25 = sapply(rtss_none, `[`, 4),
      rtss_none_q75 = sapply(rtss_none, `[`, 5),
      
      # Both vs RTSS
      both_rtss_median = sapply(both_rtss, `[`, 1),
      both_rtss_q025 = sapply(both_rtss, `[`, 2),
      both_rtss_q975 = sapply(both_rtss, `[`, 3),
      both_rtss_q25 = sapply(both_rtss, `[`, 4),
      both_rtss_q75 = sapply(both_rtss, `[`, 5),
      
      # SMC vs none
      smc_none_median = sapply(smc_none, `[`, 1),
      smc_none_q025 = sapply(smc_none, `[`, 2),
      smc_none_q975 = sapply(smc_none, `[`, 3),
      smc_none_q25 = sapply(smc_none, `[`, 4),
      smc_none_q75 = sapply(smc_none, `[`, 5),
      
      # Both vs none
      both_none_median = sapply(both_none, `[`, 1),
      both_none_q025 = sapply(both_none, `[`, 2),
      both_none_q975 = sapply(both_none, `[`, 3),
      both_none_q25 = sapply(both_none, `[`, 4),
      both_none_q75 = sapply(both_none, `[`, 5),
      
      # RTSS vs SMC
      rtss_smc_median = sapply(rtss_smc, `[`, 1),
      rtss_smc_q025 = sapply(rtss_smc, `[`, 2),
      rtss_smc_q975 = sapply(rtss_smc, `[`, 3),
      rtss_smc_q25 = sapply(rtss_smc, `[`, 4),
      rtss_smc_q75 = sapply(rtss_smc, `[`, 5),
      
      # SMC vs RTSS
      smc_rtss_median = sapply(smc_rtss, `[`, 1),
      smc_rtss_q025 = sapply(smc_rtss, `[`, 2),
      smc_rtss_q975 = sapply(smc_rtss, `[`, 3),
      smc_rtss_q25 = sapply(smc_rtss, `[`, 4),
      smc_rtss_q75 = sapply(smc_rtss, `[`, 5),
      
      # Expected efficacy
      expected_efficacy_median = sapply(expected_efficacy, `[`, 1),
      expected_efficacy_q025 = sapply(expected_efficacy, `[`, 2),
      expected_efficacy_q975 = sapply(expected_efficacy, `[`, 3),
      expected_efficacy_q25 = sapply(expected_efficacy, `[`, 4),
      expected_efficacy_q75 = sapply(expected_efficacy, `[`, 5),
      
      # Ratio pred/exp
      ratio_pred_exp_median = sapply(ratio_pred_exp, `[`, 1),
      ratio_pred_exp_q025 = sapply(ratio_pred_exp, `[`, 2),
      ratio_pred_exp_q975 = sapply(ratio_pred_exp, `[`, 3),
      ratio_pred_exp_q25 = sapply(ratio_pred_exp, `[`, 4),
      ratio_pred_exp_q75 = sapply(ratio_pred_exp, `[`, 5)
    ) %>%
    select(-both_smc, -rtss_none, -both_rtss, -smc_none, -both_none, 
           -expected_efficacy, -ratio_pred_exp)
  
  # Make long 
  inci_long <- inci_summary %>%
    select(time_value, time_value_num, time_unit, starts_with('incidence')) %>%
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
    select(time_value, time_value_num, time_unit, both_smc_median:smc_rtss_q75, expected_efficacy_median:expected_efficacy_q75) %>%
    pivot_longer(cols =  c(both_smc_median:smc_rtss_q75,expected_efficacy_median:expected_efficacy_q75),
                 names_to = 'comparison',
                 values_to = "irr")%>%
    separate(comparison, 
             into = c('comparison','statistic'),
             sep = '_(?=[^_]+$)') %>%
    mutate(
          # Clean up the comparison names for better labels
          comparison = gsub("_", " vs ", comparison),
          comparison = ifelse(comparison == 'expected vs efficacy', 'Expected both vs none', comparison)
        ) %>%
    pivot_wider(
      names_from = statistic,
      values_from = irr
    )%>% 
    mutate(metric = 'efficacy',
           arm = NA)
  
  exp_ratio_long <- inci_summary %>%
    select(time_value,time_value_num,  time_unit, contains('ratio')) %>%
    pivot_longer(cols =  contains('ratio'),
                 names_to = 'statistic',
                 values_to = "ratio_pred_exp",
                 names_prefix = "ratio_pred_exp_") %>%
    pivot_wider(
      names_from = statistic,
      values_from = ratio_pred_exp
    ) %>% 
    mutate(metric = 'ratio pred to exp', arm = NA, comparison = NA)
  
  inci_summary_all <- rbind(inci_long, 
                            irrs_long,
                            exp_ratio_long) %>%
    rename(lower_ci = q025, 
           upper_ci = q975)
  
  saveRDS(inci_summary, paste0(path, outputsfolder, '/inci_summary_wide_', agg_unit, '.rds'))
  saveRDS(inci_summary_all, paste0(path, outputsfolder, '/inci_summary_all_', agg_unit, '.rds'))
}