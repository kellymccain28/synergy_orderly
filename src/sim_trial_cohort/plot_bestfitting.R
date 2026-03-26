# Plot figures for chapter 6 best fitting 
# need to have already run the make_figures.R on the repetitions in the outputs/ folder that uses the coefs from the optimisation (outputs_fitting/)

plot_bestfitting <- function(country_to_use){
  path = 'R:/Kelly/synergy_orderly/src/'
  # Pull in repetitions of the best-fitting coefficients 
  if(country_to_use == 'Mali'){ 
    optimfolder = 'sim_trial_cohort/outputs_fitting/outputs_2026-03-23_Mali'
    folder <- 'sim_trial_cohort/outputs/outputs_2026-03-24_2' ## update this when optimisation is done and repetitions have been run
    incidence_trial <- readRDS('R:/Kelly/synergy_orderly/archive/trial_results/20260219-104643-cb128f65/monthly_incidence_trial_Mali.rds')
  } else if (country_to_use == 'BF') { 
    optimfolder = 'sim_trial_cohort/outputs_fitting/outputs_2026-03-23_BF'
    folder <- 'sim_trial_cohort/outputs/outputs_2026-03-24' ## update this when optimisation is done and repetitions have been run
    incidence_trial <- readRDS('R:/Kelly/synergy_orderly/archive/trial_results/20260219-104643-cb128f65/monthly_incidence_trial_BF.rds')
  } 
  
  inci_model <- readRDS(paste0(path, folder, '/incidence.rds'))
  infs_formatted_model <- readRDS(paste0(path, folder, '/formatted_infrecords.rds'))
  
  # Pull in trial summary results 
  tidy_results_trial <- readRDS('R:/Kelly/synergy_orderly/archive/trial_results/20260324-110905-235eebeb/surv_analysis_trial_stratified.rds') %>%
    filter(country == country_to_use)
  
  
  # Summarize the repetitions
  incidence_model_summ <- inci_model %>%
    group_by(arm, year, month, yearmonth) %>%
    # # For each simulation draw -- i was testing the below c(commented out) to see if there was a diff between the CIs i did on the central inci and the PIs from simulation -- very similar 
    # ggplot(incidence_model_summ) + geom_line(aes(x = yearmonth, y = incidence_param_median, color = 'param uncert')) +geom_line(aes(x = yearmonth, y = incidence_pred_median, color = 'pred uncert'))+facet_wrap(~arm)
    # mutate(
    #   # Simulate new cases from the predicted rate
    #   sim_cases = rpois(n(), lambda = incidence_per_1000pm * person_months / 1000),
    #   sim_incidence = (sim_cases / person_months) * 1000
    # ) %>%
    # # Then summarize predictive distribution
    # summarize(
    #   incidence_pred_median = median(sim_incidence, na.rm = TRUE),
    #   incidence_pred_lower = quantile(sim_incidence, 0.025, na.rm = TRUE),
    #   incidence_pred_upper = quantile(sim_incidence, 0.975, na.rm = TRUE),
    #   # Keep parameter uncertainty version too for comparison
    #   incidence_param_median = median(incidence_per_1000pm, na.rm = TRUE),
    #   incidence_param_lower = quantile(incidence_per_1000pm, 0.025, na.rm = TRUE),
    #   incidence_param_upper = quantile(incidence_per_1000pm, 0.975, na.rm = TRUE)
    # ) %>%
    summarize(across(c(incidence_per_1000pm, person_months, n_cases),
                     list(median = ~quantile(.x, 0.5, na.rm = TRUE),
                          lower = ~quantile(.x, 0.025, na.rm = TRUE),
                          upper = ~quantile(.x, 0.975, na.rm = TRUE))
    ) ) %>%
    # rename those variables with _median to be just the variable name 
    rename_with(.fn = \(x)sub("_median","", x)) %>%
    select(yearmonth, arm, 
           incidence_per_1000pm_model = incidence_per_1000pm,
           incidence_per_1000pm_model_lower = incidence_per_1000pm_lower,
           incidence_per_1000pm_model_upper = incidence_per_1000pm_upper,
           person_months_model = person_months, 
           n_cases_model = n_cases)
  
  incidence_trial <- incidence_trial %>%
    select(yearmonth, arm, 
           incidence_per_1000pm_trial = incidence_per_1000pm,
           person_months_trial = person_months, 
           n_cases_trial = n_cases,
           incidence_per_1000pm_trial_lower = lower_per_1000, 
           incidence_per_1000pm_trial_upper = upper_per_1000)
  
  inci_joined <- left_join(incidence_model_summ,
                           incidence_trial) %>%
    filter(arm != 'none') %>%
    mutate(yearmonth = as.Date(yearmonth))
  
  # Plot the comparison of the incidence in trial and in model output
  inc <- ggplot(inci_joined)+
    geom_point(aes(x = as.Date(yearmonth), y = incidence_per_1000pm_trial, color = 'Trial'), size= 0.9) +
    geom_point(aes(x = as.Date(yearmonth), y = incidence_per_1000pm_model, color = 'Model'), size = 0.9) +
    # geom_line(aes(x = yearmonth, y = incidence_per_1000pm_trial, color = 'Trial'), linewidth = 0.4) +
    # geom_line(aes(x = yearmonth, y = incidence_per_1000pm_model, color = 'Model'), linewidth = 0.4) +
    geom_errorbar(aes(x = as.Date(yearmonth), ymin = incidence_per_1000pm_trial_lower, ymax = incidence_per_1000pm_trial_upper, color = 'Trial'),
                  alpha  = 1, width = 15, linewidth = 0.25) +
    geom_errorbar(aes(x = as.Date(yearmonth), ymin = incidence_per_1000pm_model_lower, ymax = incidence_per_1000pm_model_upper, color = 'Model'),
                alpha  = 0.9, width = 15, linewidth = 0.25) +
    scale_color_manual(values =  c('Trial' = '#E15554', 
                                   'Model' = '#E1BC29'))+#c('#C44536','#772E25','#197278','#283D3B'))+
    scale_x_date(breaks = '3 months',
                 labels = scales::label_date_short()) +
    labs(color = NULL,#'Intervention arm',
         fill = NULL,#'Intervention arm',
         x = 'Date',
         y = 'Incidence per 1000 person-months') +
    theme_minimal(base_size = 14) + 
    facet_wrap(~factor(arm, levels = c('none','rtss','smc','both')), 
               nrow = 4,
               scales = 'free')
  inc
  # Plot efficacy comparison between model and trial 
  # (within get_cox_efficacy, I deal with the multiple simulations)
  eff_model_smc <- get_cox_efficacy(df = infs_formatted_model, 
                                    ref = 'arm_smcref',
                                    model = TRUE)
  eff_model_rtss <- get_cox_efficacy(df = infs_formatted_model, 
                                     ref = 'arm_rtssref',
                                     model = TRUE)
  model_results <- rbind(eff_model_smc, eff_model_rtss) %>%
    mutate(year = factor(year, 
                         levels = c(1, 2, 3, 'overall'),
                         labels = c("Year 1", "Year 2", "Year 3", "Overall")))
  
  eff <- ggplot(tidy_results_trial %>% 
                         filter(term %in% c("Both vs. RTSS",'RTSS vs. SMC','Both vs. SMC')))+
    geom_point(aes(x = term, y = VE, group = year, color = year, shape = 'Trial'),
               position = position_dodge(width=0.3), size = 2, alpha = 0.7) + 
    geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper, group = year, color = year, linetype = 'Trial'), 
                  position = position_dodge(width=0.3), width = 0.35, alpha = 0.7) +
    
    geom_point(data = model_results %>% 
                 filter(term %in% c("Both vs. RTSS",'RTSS vs. SMC','Both vs. SMC')),
               aes(x = term, y = VE, group = year, color = year, shape = 'Model'),
               position = position_dodge(width=0.3), size = 2) + 
    geom_errorbar(data = model_results %>% 
                    filter(term %in% c("Both vs. RTSS",'RTSS vs. SMC','Both vs. SMC')),
                  aes(x = term, ymin = VE_lower, ymax = VE_upper, group = year, color = year, linetype = 'Model'), 
                  position = position_dodge(width=0.3), width = 0.35) +
    geom_hline(aes(yintercept = 0), linetype = 2) +
    scale_y_continuous(breaks = seq(-1, 1, 0.2),
                       # limits = c(-01, 1),
                       labels = scales::percent) +
    scale_color_manual(values = c('Year 1' = '#7FB800',
                                  'Year 2' = '#573280',
                                  'Year 3' = '#CE6479',
                                  'Overall' = '#197BBD')) +
    labs(x = '',
         y = 'Efficacy',
         color = 'Trial year',
         shape = NULL, linetype = NULL) +
    theme_minimal(base_size = 14) + 
    theme(#legend.position = c(0.7, 0.8),
      legend.background = element_rect(fill = "transparent", color = NA)) 
  
  p <- plot_grid(inc, eff, nrow = 1, labels = 'AUTO', rel_widths = c(1.2,1))
  
  ggsave(paste0('R:/Kelly/synergy_orderly/src/sim_trial_cohort/thesis_plots/best_fitting_', country_to_use, '.pdf'),
                height = 6, width = 13)
}

plot_bestfitting(country_to_use = 'BF')
plot_bestfitting(country_to_use = 'Mali')
