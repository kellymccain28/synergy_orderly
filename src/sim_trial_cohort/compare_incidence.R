# Function to calculate the root mean least squared difference (2-norm) between model-predicted and observed incidence 
# loss function to minimise when fitting lag and scaler 

compare_incidence <- function(incidence_model, 
                              incidence_trial, 
                              output_dir){
  
  parameter_grid <- readRDS(paste0(output_dir, '/parameter_grid.rds'))
  pars <- parameter_grid %>% filter(sim_id == incidence_model$sim_id[1])
  # Model incidence should be from a single parameter set
  # incidence_model <- readRDS("R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs/outputs_2026-01-27_8/incidence.rds")
  # incidence_trial <- readRDS("R:/Kelly/synergy_orderly/archive/trial_results/20260127-164922-12477275/monthly_incidence_trial.rds")
  library(zoo)
  incidence_model_summ <- incidence_model %>%
    group_by(arm, year, month, yearmonth) %>%
    summarize(across(c(incidence_per_1000pm, person_months, n_cases),
                     list(median = ~quantile(.x, 0.5, na.rm = TRUE))
                     ) ) %>%
    # rename those variables with _median to be just the variable name 
    rename_with(.fn = \(x)sub("_median","", x)) %>%
    select(yearmonth, arm, 
           incidence_per_1000pm_model = incidence_per_1000pm,
           person_months_model = person_months, 
           n_cases_model = n_cases)
  
  incidence_trial <- incidence_trial %>%
    select(yearmonth, arm, 
           incidence_per_1000pm_trial = incidence_per_1000pm,
           person_months_trial = person_months, 
           n_cases_trial = n_cases)
    
  inci_joined <- left_join(incidence_model_summ,
                           incidence_trial) %>%
    filter(arm != 'none')
  message('joined model and trial incidence')
  
  # Calculate root mean squared error OVERALL -- minimise  ----
  rmse <- sqrt(mean((inci_joined$incidence_per_1000pm_trial - inci_joined$incidence_per_1000pm_model)^2))
  message('got rmse: ', rmse)
  
  # === PER-ARM RMSE ====
  rmse_by_arm <- inci_joined %>%
    group_by(arm) %>%
    summarise(
      rmse = sqrt(mean((incidence_per_1000pm_trial - incidence_per_1000pm_model)^2)),
      .groups = 'drop'
    )
  message('RMSE by arm:\n', paste(capture.output(print(rmse_by_arm)), collapse = '\n'))
  
  # === PER-ARM OVER TIME (monthly RMSE) ====
  rmse_by_arm_time <- inci_joined %>%
    group_by(arm, yearmonth) %>%
    summarise(
      rmse_monthly = sqrt(mean((incidence_per_1000pm_trial - incidence_per_1000pm_model)^2)),
      incidence_trial = mean(incidence_per_1000pm_trial),
      incidence_model = mean(incidence_per_1000pm_model),
      .groups = 'drop'
    )
  
  # Calculate the negative likelihood
  # Get model predicted rate (cases per person-month)
  model_rate <- inci_joined$n_cases_model / inci_joined$person_months_model
  
  # Expected cases given trial exposure
  expected_cases <- model_rate * inci_joined$person_months_trial
  epsilon <- 1e-6
  expected_cases <- pmax(expected_cases, epsilon)
  # inci_joined <- inci_joined %>%
  #   mutate(n_cases_model = pmax(n_cases_model, epsilon))
  
  # POISSON NEG LL 
  log_lik <- sum(dpois(inci_joined$n_cases_trial, lambda = expected_cases, log = TRUE))
  negll <- -log_lik
  
  # For Poisson with no estimated parameters (lambda is fixed from model predictions)
  # Number of parameters = 0 (the model gives expected cases directly)
  n_params <- 0
  aic <- 2 * n_params + 2 * negll  # = 2 * negll
  bic <- log(nrow(inci_joined)) * n_params + 2 * negll  # = 2 * negll
  message('got negll: ', round(negll, 2))
  message('AIC: ', round(aic, 2))
  message('BIC: ', round(bic, 2))
  
  negll_by_arm <- inci_joined %>%
    group_by(arm) %>%
    summarise({
      expected_cases_arm <- (n_cases_model / person_months_model) * person_months_trial
      expected_cases_arm <- pmax(expected_cases_arm, epsilon)
      negll_val <- -sum(dpois(n_cases_trial, lambda = expected_cases_arm, log = TRUE))
      tibble(negll = negll_val)
    }, .groups = 'drop')
  message('NegLL by arm:\n', paste(capture.output(print(negll_by_arm)), collapse = '\n'))
  
  # Per-arm Poisson metrics
  poisson_by_arm <- inci_joined %>%
    group_by(arm) %>%
    summarise({
      expected_cases_arm <- (n_cases_model / person_months_model) * person_months_trial
      expected_cases_arm <- pmax(expected_cases_arm, epsilon)
      
      log_lik_arm <- sum(dpois(n_cases_trial, 
                               lambda = expected_cases_arm, 
                               log = TRUE))
      negll_arm <- -log_lik_arm
      
      # For Poisson with fixed lambda (no estimated parameters)
      n_params_arm <- 0
      aic_arm <- 2 * negll_arm
      bic_arm <- 2 * negll_arm
      
      tibble(
        negll = negll_arm,
        aic = aic_arm,
        bic = bic_arm,
        n_obs = n(),
        rmse = sqrt(mean((incidence_per_1000pm_trial - incidence_per_1000pm_model)^2))
      )
    }, .groups = 'drop')
  
  message('\nPoisson fit by arm:')
  print(poisson_by_arm)
  
  message('got negll: ', negll)
  
  iter <- unique(incidence_model$iter)
  simid <- unique(incidence_model$sim_id)
  if(length(simid) > 1){stop("There is more than one simulation in the model incidence data frame")}
  
  # Plot the comparison of the two incidences
  incicomparison <- ggplot(inci_joined)+
    geom_line(aes(x = yearmonth, y = incidence_per_1000pm_trial, color = 'Trial'), linewidth = 0.8) +
    geom_line(aes(x = yearmonth, y = incidence_per_1000pm_model, color = 'Model'), linewidth = 0.8) +
    scale_color_manual(values =  c('Trial' = '#E15554', 
                                   'Model' = '#E1BC29'))+#c('#C44536','#772E25','#197278','#283D3B'))+
    labs(color = 'Intervention arm',
         x = 'Date',
         y = 'Incidence per 1000 person-months') +
    theme_bw(base_size = 14) + 
    facet_wrap(~factor(arm, levels = c('none','rtss','smc','both')), 
               nrow = 4,
               scales = 'free')
  
  # Plot 2: RMSE over time by arm
  rmse_time_plot <- ggplot(rmse_by_arm_time, aes(x = yearmonth, y = rmse_monthly, color = arm)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    labs(x = 'Date', 
         y = 'Monthly RMSE (incidence per 1000 pm)',
         title = paste('Model vs Trial Fit Over Time (Overall RMSE:', round(rmse, 3), ')'),
         color = 'Arm') +
    theme_bw(base_size = 12) +
    facet_wrap(~arm, scales = 'free_y', nrow = 2)
  
  # Plot 3: Scatter plot of model vs trial by arm
  scatter_plot <- ggplot(inci_joined, aes(x = incidence_per_1000pm_model, 
                                          y = incidence_per_1000pm_trial)) +
    geom_point(aes(color = arm), alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'gray50') +
    geom_smooth(method = 'lm', se = FALSE, color = 'black', linewidth = 0.5) +
    labs(x = 'Model Incidence (per 1000 pm)',
         y = 'Trial Incidence (per 1000 pm)',
         title = paste('Overall RMSE:', round(rmse, 3), 
                       ' | AIC:', round(aic, 1))) +
    theme_bw(base_size = 12) +
    facet_wrap(~arm, scales = 'free')
  
  
  # save plots 
  output_dir_plots <- paste0(output_dir, '/plots')
  if(!dir.exists(output_dir_plots)) dir.create(output_dir_plots, recursive = TRUE)
  # function to generate unique file nam e
  get_unique_filename <- function(base_path, base_name, iter, extension = '.png') {
    filename <- paste0(base_path, '/', base_name, '_', iter, extension)
    if (!file.exists(filename)) {
      return(filename)
    } else {
      counter <- 1
      repeat {
        new_filename <- paste0(base_path, '/', base_name, '_', iter, '_', counter, extension)
        if (!file.exists(new_filename)) {
          return(new_filename)
        }
        counter <- counter + 1
      }
    }
  }

  
  # Save incidence comparison plot
  incifile <- get_unique_filename(output_dir_plots, 'incidence_trial_vs_model', iter, '.png')
  ggsave(filename = incifile, incicomparison, width = 10, height = 12)
  message('Saved: ', incifile)
  
  # Save RMSE over time plot
  rmsefile <- get_unique_filename(output_dir_plots, 'rmse_over_time', iter, '.png')
  ggsave(filename = rmsefile, rmse_time_plot, width = 10, height = 8)
  message('Saved: ', rmsefile)
  
  # Save scatter plot
  scatterfile <- get_unique_filename(output_dir_plots, 'scatter_model_vs_trial', iter, '.png')
  ggsave(filename = scatterfile, scatter_plot, width = 10, height = 8)
  message('Saved: ', scatterfile)
  
  # Save RMSE metrics to CSV
  rmse_metrics <- list(
    overall = data.frame(
      iter = iter,
      simid = simid,
      rmse = rmse,
      negll = negll,
      log_lik = log_lik,
      aic = aic,
      bic = bic,
      n_obs = nrow(inci_joined)
    ),
    by_arm = rmse_by_arm,
    by_arm_time = rmse_by_arm_time,
    poisson_by_arm = poisson_by_arm
  )
  
  write.csv(rmse_metrics$overall, 
            get_unique_filename(output_dir_plots, 'metrics_overall', iter, '.csv'), 
            row.names = FALSE)
  write.csv(rmse_metrics$by_arm, 
            get_unique_filename(output_dir_plots, 'rmse_by_arm', iter, '.csv'), 
            row.names = FALSE)
  write.csv(rmse_metrics$by_arm_time, 
            get_unique_filename(output_dir_plots, 'rmse_by_arm_time', iter, '.csv'), 
            row.names = FALSE)
  write.csv(rmse_metrics$poisson_by_arm, 
            get_unique_filename(output_dir_plots, 'poisson_metrics_by_arm', iter, '.csv'), 
            row.names = FALSE)
  
  return(list(
    rmse = rmse,  
    rmse_by_arm = rmse_by_arm,
    rmse_by_arm_time = rmse_by_arm_time,
    negll = negll,  
    negll_by_arm = negll_by_arm,  # 
    log_lik = log_lik,  # 
    aic = aic,  # 
    bic = bic,  # 
    simid = simid,
    iter = iter
  ))
}

# Function to compare the hazard ratios from model and trial 
compare_hr <- function(infs_formatted_model, 
                       output_dir,
                       country_to_use){
  output_dir <- paste0(output_dir, '/plots')
  
  tidy_results_trial <- readRDS('R:/Kelly/synergy_orderly/archive/trial_results/20260324-110905-235eebeb/surv_analysis_trial_stratified.rds') %>%
    filter(country == country_to_use)
  
  iter <- unique(infs_formatted_model$iter)
  simid <- unique(infs_formatted_model$sim_id)
  if(length(simid) > 1){stop("There is more than one simulation in the formatted data frame")}
  
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
  
  # Plot efficacy
  efficacies <- ggplot(tidy_results_trial %>% 
                         filter(term %in% c("Both vs. RTSS",'RTSS vs. SMC','Both vs. SMC')))+
    geom_point(aes(x = term, y = VE, group = year, color = year, shape = 'Trial'),
               position = position_dodge(width=0.3), size = 2) + 
    geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper, group = year, color = year, linetype = 'Trial'), 
                  position = position_dodge(width=0.3), width = 0.2) +
    
    geom_point(data = model_results %>% 
                 filter(term %in% c("Both vs. RTSS",'RTSS vs. SMC','Both vs. SMC')),
               aes(x = term, y = VE, group = year, color = year, shape = 'Model'),
               position = position_dodge(width=0.3), size = 2) + 
    geom_errorbar(data = model_results %>% 
                    filter(term %in% c("Both vs. RTSS",'RTSS vs. SMC','Both vs. SMC')),
                  aes(x = term, ymin = VE_lower, ymax = VE_upper, group = year, color = year, linetype = 'Model'), 
                  position = position_dodge(width=0.3), width = 0.2) +
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
    theme_bw(base_size = 14) + 
    theme(#legend.position = c(0.7, 0.8),
      legend.background = element_rect(fill = "transparent", color = NA)) 
  
  base_filename <- paste0(output_dir, '/hr_trial_vs_model_', iter, '.png')
  
  if (!file.exists(base_filename)) {
    ggsave(filename = base_filename, efficacies)
  } else {
    # If exists, find next available number
    counter <- 1
    repeat {
      new_filename <- paste0(output_dir, '/hr_trial_vs_model_', iter, '_', counter, '.png')
      if (!file.exists(new_filename)) {
        ggsave(filename = new_filename, efficacies)
        break
      }
      counter <- counter + 1
    }
  }
  
}

# Method of moments estimator -- as an estimation of the size parameter -- may need to estimate it explicitly later 
# trial_mean <- mean(inci_joined$n_cases_trial)
# trial_var <- var(inci_joined$n_cases_trial)

# For Negative Binomial: variance = mu + mu^2/size
# Solving for size:
# size_est <- trial_mean^2 / (trial_var - trial_mean) # should be positive 
# 0.5854176

# # Look at variance-to-mean ratio by arm/month
# monthly_incidence_trial  %>%
#   group_by(arm) %>%
#   summarise(
#     mean_cases = mean(n_cases),
#     var_cases = var(n_cases),
#     ratio = var_cases / mean_cases
#   )

# best_fit <- optim(
#   par = c(param1_init, param2_init, log(1)),
#   fn = nll_negbin,  # This returns NEGATIVE log-likelihood
#   method = "L-BFGS-B",
#   lower = c(param1_lower, param2_lower, log(0.01)),
#   upper = c(param1_upper, param2_upper, log(100))
# )
# 
# # Best parameters are at:
# best_params <- best_fit$par[1:2]
# best_dispersion <- exp(best_fit$par[3])
# 
# # The MINIMUM negative log-likelihood is:
# cat("Minimum -LL:", best_fit$value, "\n")
# 
# # Check if optimization succeeded
# if(best_fit$convergence == 0) {
#   cat("Optimization converged successfully!\n")
# } else {
#   cat("Warning: convergence code", best_fit$convergence, "\n")
# }
# 
# # Look at final negative log-likelihood
# cat("Final -LL:", best_fit$value, "\n")
# 
# # Run model with best parameters
# final_expected <- run_model_with_best_params(best_params)
# 
# # Visual check
# plot(incidence_data$trial_cases, final_expected, 
#      xlab = "Observed", ylab = "Predicted",
#      main = "Model Fit")
# abline(0, 1, col = "red")
# 
# # Goodness of fit
# cor(incidence_data$trial_cases, final_expected)