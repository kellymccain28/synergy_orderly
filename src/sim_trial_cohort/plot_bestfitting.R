# Plot figures for chapter 6 best fitting 
# need to have already run the make_figures.R on the repetitions in the outputs/ folder that uses the coefs from the optimisation (outputs_fitting/)

# rmse_values <- sapply(1:200, function(ii) {
#   val <- optim_checkpoint$history[[ii]]$rmse
#   if(is.null(val)) Inf else val
# })
# 
# # Get indices ordered by RMSE
# ordered_indices <- order(rmse_values)
# 
# # Get the index of the second lowest RMSE
# second_lowest_index <- ordered_indices[2]
# 
# # Pull out the coefficients for that index
# coefs_second_lowest <- optim_checkpoint$history[[second_lowest_index]]$coefs

plot_bestfitting <- function(country_to_use,
                             num_arm_fit){
  path = 'R:/Kelly/synergy_orderly/src/'
  # Pull in repetitions of the best-fitting coefficients 
  if(country_to_use == 'Mali'){ 
    # optimfolder = 'sim_trial_cohort/outputs_fitting/outputs_2026-03-23_Mali'
    if(num_arm_fit == 3){
      folder <- 'sim_trial_cohort/outputs/outputs_2026-03-24_2' ## this is for best-fitting values when fitting to 3 arms (03-23 fitting)
    }  else if(num_arm_fit == 2){
      # folder <- 'sim_trial_cohort/outputs/outputs_2026-03-26_2' ## this is for best-fitting values when fitting to 2 arms (03-25_2 fitting - pre-finished optim)
      folder <- 'sim_trial_cohort/outputs/outputs_2026-03-30_2' ## this is for best-fitting values when fitting to 2 arms (03-25_2 fitting)
    }
    incidence_trial <- readRDS('R:/Kelly/synergy_orderly/archive/trial_results/20260219-104643-cb128f65/monthly_incidence_trial_Mali.rds')
  } else if (country_to_use == 'BF') { 
    # optimfolder = 'sim_trial_cohort/outputs_fitting/outputs_2026-03-23_BF'
    if(num_arm_fit == 3){
      folder <- 'sim_trial_cohort/outputs/outputs_2026-03-24' ## this is for best-fitting values when fitting to 3 arms (03-23 fitting)
    } else if(num_arm_fit == 2){
      # folder <- 'sim_trial_cohort/outputs/outputs_2026-03-26' ## this is for best-fitting values when fitting to 2 arms (03-25_2 fitting - pre-finished optim)
      folder <- 'sim_trial_cohort/outputs/outputs_2026-03-30' ## this is for best-fitting values when fitting to 2 arms (03-25_2 fitting)
    }
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
    geom_line(aes(x = yearmonth, y = incidence_per_1000pm_trial, color = 'Trial'), linewidth = 0.4) +
    geom_line(aes(x = yearmonth, y = incidence_per_1000pm_model, color = 'Model'), linewidth = 0.4) +
    geom_errorbar(aes(x = as.Date(yearmonth), ymin = incidence_per_1000pm_trial_lower, ymax = incidence_per_1000pm_trial_upper, color = 'Trial'),
                  alpha  = 1, width = 15, linewidth = 0.25) +
    geom_errorbar(aes(x = as.Date(yearmonth), ymin = incidence_per_1000pm_model_lower, ymax = incidence_per_1000pm_model_upper, color = 'Model'),
                  alpha  = 0.9, width = 15, linewidth = 0.25) +
    scale_color_manual(values =  c('Trial' = '#E1BC29',      
                                   'Model' = '#E15554'))+#c('#C44536','#772E25','#197278','#283D3B'))+
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
         y = 'Relative efficacy',
         color = 'Trial year',
         shape = NULL, linetype = NULL) +
    theme_minimal(base_size = 14) + 
    theme(#legend.position = c(0.7, 0.8),
      legend.background = element_rect(fill = "transparent", color = NA)) 
  
  p <- plot_grid(inc, eff, nrow = 1, labels = 'AUTO', rel_widths = c(1.2,1))
  
  ggsave(paste0('R:/Kelly/synergy_orderly/src/sim_trial_cohort/thesis_plots/best_fitting_', country_to_use, '_', num_arm_fit,'armfit.pdf'),
         height = 6, width = 13)
}

plot_bestfitting(country_to_use = 'BF',
                 num_arm_fit = 3)
plot_bestfitting(country_to_use = 'Mali',
                 num_arm_fit = 3)

plot_bestfitting(country_to_use = 'BF',
                 num_arm_fit = 2)
plot_bestfitting(country_to_use = 'Mali',
                 num_arm_fit = 2)


#  Plot splines
pbitemali3 <- readRDS('R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs_fitting/outputs_2026-03-23_Mali/best_so_far.rds')$params_row$p_bite[[1]] 
pbitebf3 <- readRDS('R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs_fitting/outputs_2026-03-23_BF/best_so_far.rds')$params_row$p_bite[[1]]
pbite3 <- data.frame(country = c(rep('Mali',length(pbitemali3)), rep('BF',length(pbitebf3))),
                     pbite = c(pbitemali3, pbitebf3),
                     date = seq(as.Date('2017-04-01')-50, as.Date('2020-03-31')),
                     armfit = '3-arm fit') 
pbitemali2 <- readRDS('R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs_fitting/outputs_2026-03-25_Mali_2/best_so_far.rds')$params_row$p_bite[[1]]
pbitebf2 <- readRDS('R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs_fitting/outputs_2026-03-25_BF_2/best_so_far.rds')$params_row$p_bite[[1]]
pbite2 <- data.frame(country = c(rep('Mali',length(pbitemali2)), rep('BF',length(pbitebf2))),
                     pbite = c(pbitemali2, pbitebf2),
                     date = seq(as.Date('2017-04-01')-50, as.Date('2020-03-31')),
                     armfit = '2-arm fit')

pbite23 <- bind_rows(pbite2, pbite3)

ggplot(pbite3) + 
  geom_point(aes(x = date, y = pbite, group = country), color = '#E15554') + 
  facet_wrap(~country, nrow = 2) + 
  labs(x = 'Date',
       y = 'Best-fitting spline for the daily probability of an infectious bite') +
  theme_bw()
ggsave(paste0('R:/Kelly/synergy_orderly/src/sim_trial_cohort/thesis_plots/best_fitting_3arm_splines.pdf'),
       height = 6, width = 8)

ggplot(pbite23) + 
  geom_line(aes(x = date, y = pbite, group = armfit, color = armfit), alpha = 0.8, linewidth = 1) + 
  facet_wrap(~country, nrow = 2) + 
  scale_color_manual(values =  c('2-arm fit' = 'darkred', 
                                 '3-arm fit' = '#E15554')) + 
  labs(x = 'Date', color = NULL,
       y = 'Best-fitting spline for the daily probability of an infectious bite') +
  theme_bw()
ggsave(paste0('R:/Kelly/synergy_orderly/src/sim_trial_cohort/thesis_plots/best_fitting_2and3armfit_splines.pdf'),
       height = 6, width = 8)



# Below is code to plot the 3-arm fit in grey with the 2-arm fit in color on top (modified from function above)

path = 'R:/Kelly/synergy_orderly/src/'
outputs_folders <-  c(#'BF3syn' = 'outputs_2026-03-25_6',
                      # 'Mali3syn' = 'outputs_2026-03-25_5',
                      # 'BF2syn' = 'outputs_2026-03-30_7',
                      # 'Mali2syn' = 'outputs_2026-03-30_8',
                      
                      'BF3best' = 'outputs_2026-03-24',
                      'Mali3best' = 'outputs_2026-03-24_2',
                      'BF2best' = 'outputs_2026-03-30',
                      'Mali2best' = 'outputs_2026-03-30_2')
# Pull in repetitions of the best-fitting coefficients 
# Define a function to read both files from a folder
read_model_outputs <- function(folder_name, folder_path, base_path = "R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs/") {
  incidence <- readRDS(file.path(base_path, folder_path, 'incidence.rds'))
  formatted <- readRDS(file.path(base_path, folder_path, 'formatted_infrecords.rds'))
  return(list(incidence = incidence, formatted = formatted))
}

# Use lapply to read all folders
results <- lapply(names(outputs_folders), function(name) {
  read_model_outputs(name, outputs_folders[name])
})
names(results) <- names(outputs_folders)

# Extract and bind incidence
incidence_combined <- map_dfr(results, ~ .x$incidence, .id = "scenario") %>%
  mutate(
  country = case_when(
    grepl("^BF", scenario) ~ "BF",
    grepl("^Mali", scenario) ~ "Mali"
  ),
  armfit = case_when(
    grepl("3", scenario) ~ "3-arm fit",
    grepl("2", scenario) ~ "2-arm fit"
  ),
  type = case_when(
    grepl("syn$", scenario) ~ "syntest",
    grepl("best$", scenario) ~ "bestreps"
  )
)

# Extract and bind formatted
formatted_combined <- map_dfr(results, ~ .x$formatted, .id = "scenario") %>%
  mutate(
    country = case_when(
      grepl("^BF", scenario) ~ "BF",
      grepl("^Mali", scenario) ~ "Mali"
    ),
    armfit = case_when(
      grepl("3", scenario) ~ "3-arm fit",
      grepl("2", scenario) ~ "2-arm fit"
    ),
    type = case_when(
      grepl("syn$", scenario) ~ "syntest",
      grepl("best$", scenario) ~ "bestreps"
    )
  )

# Pull in trial summary results 
tidy_results_trial <- readRDS('R:/Kelly/synergy_orderly/archive/trial_results/20260324-110905-235eebeb/surv_analysis_trial_stratified.rds') 

# trial incidence
incidence_trialBF <- readRDS('R:/Kelly/synergy_orderly/archive/trial_results/20260219-104643-cb128f65/monthly_incidence_trial_BF.rds') %>%
  mutate(country = 'BF')
incidence_trialMali<- readRDS('R:/Kelly/synergy_orderly/archive/trial_results/20260219-104643-cb128f65/monthly_incidence_trial_Mali.rds') %>%
  mutate(country = 'Mali')
incidence_trial <- bind_rows(incidence_trialBF, incidence_trialMali)

# Summarize the repetitions
incidence_model_summ <- incidence_combined %>%
  group_by(arm, year, month, yearmonth, scenario, type, country, armfit) %>%
  summarize(across(c(incidence_per_1000pm, person_months, n_cases),
                   list(median = ~quantile(.x, 0.5, na.rm = TRUE),
                        lower = ~quantile(.x, 0.025, na.rm = TRUE),
                        upper = ~quantile(.x, 0.975, na.rm = TRUE))
  ) ) %>%
  # rename those variables with _median to be just the variable name 
  rename_with(.fn = \(x)sub("_median","", x)) %>%
  select(scenario, yearmonth, arm, type, armfit, country,
         incidence_per_1000pm_model = incidence_per_1000pm,
         incidence_per_1000pm_model_lower = incidence_per_1000pm_lower,
         incidence_per_1000pm_model_upper = incidence_per_1000pm_upper,
         person_months_model = person_months, 
         n_cases_model = n_cases)

incidence_trial <- incidence_trial %>%
  select(country, yearmonth, arm,  
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
  # Trial data
  geom_point(aes(x = as.Date(yearmonth), y = incidence_per_1000pm_trial, color = 'Trial'), size= 0.9) +
  geom_line(aes(x = yearmonth, y = incidence_per_1000pm_trial, color = 'Trial'), linewidth = 0.4, alpha = 0.8) +
  geom_errorbar(aes(x = as.Date(yearmonth), ymin = incidence_per_1000pm_trial_lower, ymax = incidence_per_1000pm_trial_upper, color = 'Trial'),
                alpha  = 1, width = 15, linewidth = 0.25) +
  geom_ribbon(aes(x = as.Date(yearmonth), ymin = incidence_per_1000pm_trial_lower, ymax = incidence_per_1000pm_trial_upper, fill = 'Trial'),
                alpha  = 0.2, width = 15, linewidth = 0.25) +
  # Model - 3-arm
  geom_point(data = inci_joined %>% filter(armfit == '3-arm fit'),
             aes(x = as.Date(yearmonth), y = incidence_per_1000pm_model, color = '3-arm fit', group = armfit), size = 0.9) +
  geom_line(data = inci_joined %>% filter(armfit == '3-arm fit'),
            aes(x = yearmonth, y = incidence_per_1000pm_model, color = '3-arm fit', group = scenario), linewidth = 0.4, alpha = 0.8) +
  geom_errorbar(data = inci_joined %>% filter(armfit == '3-arm fit'),
                aes(x = as.Date(yearmonth), ymin = incidence_per_1000pm_model_lower, ymax = incidence_per_1000pm_model_upper, 
                    color = '3-arm fit', group = scenario),
                alpha  = 0.9, width = 15, linewidth = 0.25) +
  # Model - 2-arm
  geom_point(data = inci_joined %>% filter(armfit == '2-arm fit'),
             aes(x = as.Date(yearmonth), y = incidence_per_1000pm_model, color = '2-arm fit', group = armfit), size = 0.9) +
  geom_line(data = inci_joined %>% filter(armfit == '2-arm fit'),
            aes(x = yearmonth, y = incidence_per_1000pm_model, color = '2-arm fit', group = scenario), linewidth = 0.4, alpha = 0.8) +
  geom_errorbar(data = inci_joined %>% filter(armfit == '2-arm fit'),
                aes(x = as.Date(yearmonth), ymin = incidence_per_1000pm_model_lower, ymax = incidence_per_1000pm_model_upper, 
                    color = '2-arm fit', group = scenario),
                alpha  = 0.9, width = 15, linewidth = 0.25) +
  scale_color_manual(values =  c('Trial' = '#E1BC29',      
                                 '2-arm fit' = 'darkred', 
                                 '3-arm fit' = '#E15554'))+#c('#C44536','#772E25','#197278','#283D3B'))+
  scale_fill_manual(values =  c('Trial' = '#E1BC29',      
                                 '2-arm fit' = 'darkred', 
                                 '3-arm fit' = '#E15554'))+#c('#C44536','#772E25','#197278','#283D3B'))+
  scale_x_date(breaks = '3 months',
               labels = scales::label_date_short()) +
  labs(color = NULL,#'Intervention arm',
       fill = NULL,#'Intervention arm',
       x = 'Date',
       y = 'Incidence per 1000 person-months') +
  theme_bw(base_size = 14) + 
  facet_grid(rows = vars(factor(arm, levels = c('rtss','smc','both'))), 
             cols = vars(country),
             scales = 'free') +   guides(fill = "none") 
  # facet_wrap(~factor(arm, levels = c('none','rtss','smc','both')), 
  #            nrow = 4,
  #            scales = 'free')
inc



# Plot efficacy comparison between model and trial 
trial_overall <- tidy_results_trial %>%
  filter(term %in% c("Both vs. RTSS", 'RTSS vs. SMC', 'Both vs. SMC'),
         year == 'Overall') %>%
  mutate(
    armfit = "Trial",
    source_type = "Trial"
  )

# (within get_cox_efficacy, I deal with the multiple simulations)
scenarios_list <- split(formatted_combined, formatted_combined$scenario)

model_results_list <- lapply(names(scenarios_list), function(scen) {
  df <- scenarios_list[[scen]]
  
  eff_smc <- get_cox_efficacy(df = df, ref = 'arm_smcref', model = TRUE)
  eff_rtss <- get_cox_efficacy(df = df, ref = 'arm_rtssref', model = TRUE)
  
  bind_rows(eff_smc, eff_rtss) %>%
    mutate(
      scenario = scen,
      country = unique(df$country),
      armfit = unique(df$armfit),
      type = unique(df$type)
    )
})
model_results <- bind_rows(model_results_list) %>%
  mutate(
    year = factor(year, 
                  levels = c(1, 2, 3, 'overall'),
                  labels = c("Year 1", "Year 2", "Year 3", "Overall")),
    source_type = "Model"
  ) %>%
  filter(term %in% c("Both vs. RTSS", 'RTSS vs. SMC', 'Both vs. SMC') &
           year == 'Overall')

# Combine both datasets
combined_overall <- bind_rows(trial_overall, model_results) %>%
  mutate(
    # Create a dodge group combining term and armfit
    dodge_group = interaction(term, armfit)
  )

# Now create the plot with both fit types
eff <- ggplot(combined_overall, aes(x = term, y = VE)) +
  # Points with dodging
  geom_point(aes(color = armfit, group = dodge_group),
             position = position_dodge(width = 0.3), 
             size = 2) +
  # Error bars with dodging
  geom_errorbar(aes(ymin = VE_lower, ymax = VE_upper, 
                    color = armfit, group = dodge_group),
                position = position_dodge(width = 0.3),
                width = 0.35) +
  # Horizontal line at 0
  geom_hline(aes(yintercept = 0), linetype = 2) +
  # Scale customization
  scale_y_continuous(breaks = seq(-1, 1, 0.2),
                     labels = scales::percent) +
  scale_color_manual(name = NULL,
                     values = c('Trial' = '#E1BC29',      
                                '2-arm fit' = 'darkred', 
                                '3-arm fit' = '#E15554')) + 
  labs(x = '',
       y = 'Relative efficacy') +
  theme_bw(base_size = 14) +
  theme(legend.background = element_rect(fill = "transparent", color = NA),
        legend.position = "bottom") +
  # Facet by country
  facet_wrap(~country, ncol = 2)

eff


ggsave(paste0('R:/Kelly/synergy_orderly/src/sim_trial_cohort/thesis_plots/best_fitting_incidence_2and3armfits.pdf'),
       inc, height = 8, width = 12)
ggsave(paste0('R:/Kelly/synergy_orderly/src/sim_trial_cohort/thesis_plots/best_fitting_efficacy_2and3armfits.pdf'),
       eff, height = 5, width = 8)
