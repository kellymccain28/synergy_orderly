# Function to calculate the men least squared difference (2-norm) between model-predicted and observed incidence 
# loss function to minimise when fitting lag and scaler 
compare_incidence <- function(incidence_model, 
                              incidence_trial){
  
  # incidence_model <- readRDS("R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs/outputs_2026-01-27_8/incidence.rds")
  # incidence_trial <- readRDS("R:/Kelly/synergy_orderly/archive/trial_results/20260127-164922-12477275/monthly_incidence_trial.rds")
  
  incidence_model <- incidence_model %>%
    select(yearmonth, arm, 
           incidence_per_1000pm_model = incidence_per_1000pm,
           person_months_model = person_months, 
           n_cases_model = n_cases)
  
  incidence_trial <- incidence_trial %>%
    select(yearmonth, arm, 
           incidence_per_1000pm_trial = incidence_per_1000pm,
           person_months_trial = person_months, 
           n_cases_trial = n_cases)
    
  inci_joined <- left_join(incidence_model,
                           incidence_trial) %>%
    filter(arm != 'none')
  
  # Calculate root mean squared error 
  rmse <- sqrt(mean((inci_joined$incidence_per_1000pm_trial - inci_joined$incidence_per_1000pm_model)^2))
  
  # Calculate the negative binomial negative likelihood
  # Get model predicted rate (cases per person-month)
  # model_rate <- inci_joined$n_cases_model / inci_joined$person_months_model
  
  # Expected cases given trial exposure
  # expected_cases <- model_rate * inci_joined$person_months_trial
  epsilon <- 1e-6
  # expected_cases <- pmax(expected_cases, epsilon)
  inci_joined <- inci_joined %>%
    mutate(n_cases_model = pmax(n_cases_model, epsilon))
  
  size = 0.5854176
  negll <- -sum(dnbinom(inci_joined$n_cases_trial, 
                        mu = inci_joined$n_cases_model,#expected_cases, 
                        size = size, log = TRUE))
  
  
  return(list(
    rmse = rmse,
    negll = negll
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