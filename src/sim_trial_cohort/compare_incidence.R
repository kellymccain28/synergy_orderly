# Function to calculate the men least squared difference (2-norm) between model-predicted and observed incidence 
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
  # incidence_model <- incidence_model 
  incidence_model_summ <- incidence_model %>%
    group_by(arm, year, month, yearmonth) %>%
    summarize(across(c(incidence_per_1000pm, person_months, n_cases),
                     list(#lower = ~quantile(.x, 0.025, na.rm = TRUE),
                          median = ~quantile(.x, 0.5, na.rm = TRUE)#,
                          #upper = ~quantile(.x, 0.975, na.rm = TRUE)
                          )#,
                     # .names = "{.col}_{.fn}"
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
  
  # Calculate root mean squared error -- minimise 
  rmse <- sqrt(mean((inci_joined$incidence_per_1000pm_trial - inci_joined$incidence_per_1000pm_model)^2))
  message('got rmse: ', rmse)
  # Calculate the negative binomial negative likelihood
  # Get model predicted rate (cases per person-month)
  model_rate <- inci_joined$n_cases_model / inci_joined$person_months_model
  
  # Expected cases given trial exposure
  expected_cases <- model_rate * inci_joined$person_months_trial
  epsilon <- 1e-6
  expected_cases <- pmax(expected_cases, epsilon)
  # inci_joined <- inci_joined %>%
  #   mutate(n_cases_model = pmax(n_cases_model, epsilon))
  
  size = 0.5854176
  negll <- -sum(dnbinom(inci_joined$n_cases_trial, 
                        mu = expected_cases, #inci_joined$n_cases_model,#
                        size = size, log = TRUE))
  message('got negll: ', negll)
  
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
         y = 'Incidence per 1000 person-months',
         caption = paste0('lag: ', pars$lag_p_bite, 'yr 1 scaler: ', pars$p_bite_scaler_1,'\nyr 2 scaler: ', pars$p_bite_scaler_2, 'yr3 scaler: ',pars$p_bite_scaler_3)) +
    theme_bw(base_size = 14) + 
    facet_wrap(~factor(arm, levels = c('none','rtss','smc','both')), 
               nrow = 4)
  ggsave(filename = paste0(output_dir, '/incidence_trial_vs_model_', simid, '.png'), incicomparison)
  
  return(list(
    rmse = rmse,
    negll = negll,
    simid = simid
  ))
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