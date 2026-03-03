library(splines)
prob_bitebf <- unlist(readRDS("R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs/outputs_2026-02-19_2/parameter_grid.rds")[42,]$p_bite)
prob_biteml <- unlist(readRDS("R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs/outputs_2026-03-02/parameter_grid.rds")[1,]$p_bite)
  
# 1 Create spline from existing curve 
# BF
# pbite <- prob_bitebf %>% filter(date >= (as.Date('2017-04-01')-50) & date <= '2020-03-31') %>% pull(prob_infectious_bite)
tt <- seq(0, 365*3 + 50, 1)  # your time vector
selected_df <- 15

# Create basis matrices ONCE and save them
ns_basis <- ns(tt, df = selected_df)
X <- cbind(1, predict(ns_basis, tt))  # design matrix

# # Get starting coefficients from your curve
# starting_coefs <- coef(lm(pbite ~ ns_basis))
# Transform starting probabilities to logit scale
# Add small epsilon to avoid 0 or 1
epsilon <- 1e-6
p_safe <- pmax(pmin(prob_bitebf, 1 - epsilon), epsilon)
logit_p <- log(p_safe / (1 - p_safe))

# Fit spline on logit scale 
fit_bf <- lm(logit_p ~ ns_basis)

# Get starting coefficients on logit scale
starting_coefs <- coef(fit_bf)

# Get fitted values (on prob scale)
fitted_logit <- predict(fit_bf)
# Transform back from logit scale
fitted_logit2 <- as.numeric(X %*% starting_coefs) # this is the same as the predict(fit_ml), just doesn't use model obj
fitted_prob_bf <- plogis(fitted_logit) # # inverse logit so now in prob scale

# RMSE
rmse_bf <- sqrt(mean((prob_bitebf - fitted_prob_bf)^2))


plot(tt, prob_bitebf, type = "l", col = "blue", lwd = 1,
     main = paste("BF - Spline Fit (df =", selected_df, ")",
                  "\nRMSE =", round(rmse_bf, 6)),
     xlab = "Time (days)", ylab = "P(bite)")
lines(tt, fitted_prob_bf, col = "red", lwd = 2)
legend("topright", legend = c("Original", "Spline"), 
       col = c("blue", "red"), lty = 1, lwd = c(1, 2))

par(mfrow = c(1, 1))

# Save everything to use later
saveRDS(list(
  X = X,
  ns_basis = ns_basis,
  tt = tt,
  starting_coefs = starting_coefs,
  fitted_prob = fitted_prob_bf,
  selected_df = selected_df
), file = "R:/Kelly/synergy_orderly/src/sim_trial_cohort/spline_setupBF.rds")


# #####---------------------------------------------------------------------------------------------------
# Mali ---------------------------------------------------------------------------------------------------
# pbite <- prob_biteml %>% filter(date >= (as.Date('2017-04-01')-50) & date <= '2020-03-31')%>% pull(prob_infectious_bite)
tt <- seq(0, 365*3 + 50, 1)  # your time vector
selected_df <- 15

# Create basis matrices ONCE and save them
ns_basis <- ns(tt, df = selected_df)
X <- cbind(1, predict(ns_basis, tt))  # design matrix

# # Get starting coefficients from your curve
# starting_coefs <- coef(lm(pbite ~ ns_basis))
# Transform starting probabilities to logit scale
# Add small epsilon to avoid 0 or 1
epsilon <- 1e-6
p_safe <- pmax(pmin(prob_biteml, 1 - epsilon), epsilon)
logit_p <- log(p_safe / (1 - p_safe))

# Fit spline on logit scale 
fit_ml <- lm(logit_p ~ ns_basis)

# Get starting coefficients on logit scale
starting_coefs <- coef(fit_ml)

# Get fitted values (on prob scale)
fitted_logit <- predict(fit_ml)
# Transform back from logit scale
fitted_logit2 <- as.numeric(X %*% starting_coefs) # this is the same as the predict(fit_ml), just doesn't use model obj
fitted_prob_ml <- plogis(fitted_logit) # # inverse logit so now in prob scale

# RMSE
rmse_ml <- sqrt(mean((prob_biteml - fitted_prob_ml)^2))


plot(tt, prob_biteml, type = "l", col = "blue", lwd = 1,
     main = paste("Mali - Spline Fit (df =", selected_df, ")",
                  "\nRMSE =", round(rmse_ml, 6)),
     xlab = "Time (days)", ylab = "P(bite)")
lines(tt, fitted_prob_ml, col = "red", lwd = 2)
legend("topright", legend = c("Original", "Spline"), 
       col = c("blue", "red"), lty = 1, lwd = c(1, 2))

par(mfrow = c(1, 1))

# Save everything to use later
saveRDS(list(
  X = X,
  ns_basis = ns_basis,
  tt = tt,
  starting_coefs = starting_coefs,
  fitted_prob = fitted_prob_ml,
  selected_df = selected_df
), file = "R:/Kelly/synergy_orderly/src/sim_trial_cohort/spline_setupMali.rds")


# Load your saved setup
# setupBF <- readRDS("spline_setupBF.rds")
# setupMali <- readRDS("spline_setupMali.rds")
# XBF <- setupBF$X
# XMali <- setupMali$X
# starting_coefsBF <- setupBF$starting_coefs
# starting_coefsMali <- setupMali$starting_coefs
# 
# # Function to get probability curve (fast)
# get_prob <- function(coefs) {
#   prob <- as.numeric(X %*% coefs)
#   pmax(pmin(prob, 0.5), 0)  # bound
# }
# 
# # Wrapper for your expensive model
# run_expensive_model <- function(bite_prob) {
#   # Your existing model code here
#   # Returns incidence vector
# }
# 
# # Objective function (what optim will call)
# objective <- function(coefs) {
#   # Get probability curve
#   bite_prob <- get_prob(coefs)
#   
#   # Run expensive model
#   model_output <- run_expensive_model(bite_prob)
#   
#   # Calculate RMSE
#   rmse <- sqrt(mean((model_output - trial_data$incidence)^2))
#   
#   return(rmse)
# }
# 
# 
# # Option A: Use nloptr with limited iterations
# library(nloptr)
# 
# # Set up options for derivative-free optimization
# opts <- list(
#   algorithm = "NLOPT_LN_BOBYQA",  # Derivative-free
#   maxeval = 50,  # Only 50 evaluations! (50 * 10 min = 500 min = ~8 hours)
#   ftol_rel = 1e-3  # Stop when relative improvement is small
# )
# 
# # Run optimization
# result <- nloptr(
#   x0 = starting_coefs,
#   eval_f = objective,
#   lb = rep(-Inf, length(starting_coefs)),
#   ub = rep(Inf, length(starting_coefs)),
#   opts = opts
# )
# 
# # Option B: Use Bayesian optimization (even more efficient)
# library(mlrMBO)
# 
# # Create objective function for mlrMBO
# obj_fun <- makeSingleObjectiveFunction(
#   name = "spline_optim",
#   fn = function(x) {
#     objective(as.numeric(x))
#   },
#   par.set = makeNumericParamSet(
#     len = length(starting_coefs),
#     lower = starting_coefs - 2,  # Reasonable bounds
#     upper = starting_coefs + 2,
#     vector = TRUE
#   )
# )
# 
# # Run Bayesian optimization (usually needs fewer evaluations)
# design <- generateDesign(10, getParamSet(obj_fun), fun = lhs::randomLHS)
# control <- makeMBOControl()
# control <- setMBOControlTermination(control, max.evals = 40)
# 
# run <- mbo(obj_fun, design = design, control = control)
# 
# 
# # # Optimisation to adjust coefs
# # # Function that takes spline coefficients and returns probability curve
# # get_prob_from_coefs <- function(coefs, basis_matrix_with_intercept) {
# #   # coefs[1] is intercept, coefs[2:21] are spline coefficients
# #   prob <- basis_matrix_with_intercept %*% coefs
# #   
# #   # Ensure probabilities are in valid range (0 to 1)
# #   prob <- pmax(prob, 0)  # no negative probabilities
# #   prob <- pmin(prob, 0.5)  # reasonable upper bound
# #   
# #   return(as.numeric(prob))
# # }
# # 
# # # RMSE function for optimization
# rmse_function <- function(coefs,
#                           incidence_model,
#                           incidence_trial,
#                           output_dir) {
#   # Get probability curve from these coefficients
#   bite_prob <- get_prob_from_coefs(coefs)
# 
#   parameter_grid <- readRDS(paste0(output_dir, '/parameter_grid.rds'))
#   pars <- parameter_grid %>% filter(sim_id == incidence_model$sim_id[1])
#   # Model incidence should be from a single parameter set
#   # incidence_model <- readRDS("R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs/outputs_2026-01-27_8/incidence.rds")
#   # incidence_trial <- readRDS("R:/Kelly/synergy_orderly/archive/trial_results/20260127-164922-12477275/monthly_incidence_trial.rds")
#   library(zoo)
#   incidence_model_summ <- incidence_model %>%
#     group_by(arm, year, month, yearmonth) %>%
#     summarize(across(c(incidence_per_1000pm, person_months, n_cases),
#                      list(median = ~quantile(.x, 0.5, na.rm = TRUE))
#     ) ) %>%
#     # rename those variables with _median to be just the variable name
#     rename_with(.fn = \(x)sub("_median","", x)) %>%
#     select(yearmonth, arm,
#            incidence_per_1000pm_model = incidence_per_1000pm,
#            person_months_model = person_months,
#            n_cases_model = n_cases)
# 
#   incidence_trial <- incidence_trial %>%
#     select(yearmonth, arm,
#            incidence_per_1000pm_trial = incidence_per_1000pm,
#            person_months_trial = person_months,
#            n_cases_trial = n_cases)
# 
#   inci_joined <- left_join(incidence_model_summ,
#                            incidence_trial) %>%
#     filter(arm != 'none')
#   message('joined model and trial incidence')
# 
#   # Calculate root mean squared error -- minimise
#   rmse <- sqrt(mean((inci_joined$incidence_per_1000pm_trial - inci_joined$incidence_per_1000pm_model)^2))
#   message('got rmse: ', rmse)
# 
#   simid <- unique(incidence_model$sim_id)
#   if(length(simid) > 1){stop("There is more than one simulation in the model incidence data frame")}
# 
#   return(rmse)
# }
# # 
# # # Start from your existing coefficients
# # initial_coefs <- starting_coefs
# # 
# # # Optimize
# # result <- optim(par = initial_coefs, 
# #                 fn = rmse_function,
# #                 method = "BFGS",  # or "L-BFGS-B" if you want bounds
# #                 control = list(maxit = 1000, trace = 1))
# # 
# # # Get your fitted curve
# # fitted_coefs <- result$par
# # fitted_pbite <- get_prob_from_coefs(fitted_coefs)