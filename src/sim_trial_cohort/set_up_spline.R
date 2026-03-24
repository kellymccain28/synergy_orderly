library(splines)
prob_bitebf <- unlist(readRDS("R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs/outputs_2026-02-19_2/parameter_grid.rds")[42,]$p_bite)
prob_biteml <- unlist(readRDS("R:/Kelly/synergy_orderly/src/sim_trial_cohort/outputs/outputs_2026-03-02/parameter_grid.rds")[1,]$p_bite)
  
# 1 Create spline from existing curve 
# BF
# pbite <- prob_bitebf %>% filter(date >= (as.Date('2017-04-01')-50) & date <= '2020-03-31') %>% pull(prob_infectious_bite)
tt <- seq(0, 365*3 + 50, 1)  # time vector
selected_df <- 30

# Create basis matrix
ns_basis <- ns(tt, df = selected_df)

# Design matrix with intercept
X <- cbind(1, ns_basis)#predict(ns_basis, tt))  # design matrix

# # Get starting coefficients from the curve
# starting_coefs <- coef(lm(pbite ~ ns_basis))
# Transform starting probabilities to logit scale
# Add small epsilon to avoid 0 or 1
epsilon <- 1e-6
p_safe <- pmax(pmin(prob_bitebf, 1 - epsilon), epsilon)
logit_p <- log(p_safe / (1 - p_safe))

# Fit spline on logit scale 
fit_bf <- lm(logit_p ~ ns_basis)

# Get fitted values (on prob scale)
fitted_logit <- predict(fit_bf) # this is the same as X %*% coefficients 

# Manual for verification:
# Get starting coefficients on logit scale
starting_coefs <- coef(fit_bf)
# Transform back from logit scale
fitted_logit2 <- X %*% starting_coefs # this is the same as the predict(fit_ml), just doesn't use model obj

# Transform to probability scale 
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
), file = paste0("R:/Kelly/synergy_orderly/src/sim_trial_cohort/spline_setupBF_", selected_df, ".rds"))


# #####---------------------------------------------------------------------------------------------------
# Mali ---------------------------------------------------------------------------------------------------
# pbite <- prob_biteml %>% filter(date >= (as.Date('2017-04-01')-50) & date <= '2020-03-31')%>% pull(prob_infectious_bite)
tt <- seq(0, 365*3 + 50, 1)  # your time vector
selected_df <- 30

# Create basis matrices ONCE and save them
ns_basis <- ns(tt, df = selected_df)
X <- cbind(1, ns_basis)  # design matrix

# # Get starting coefficients from your curve
# starting_coefs <- coef(lm(pbite ~ ns_basis))
# Transform starting probabilities to logit scale
# Add small epsilon to avoid 0 or 1
epsilon <- 1e-6
p_safe <- pmax(pmin(prob_biteml, 1 - epsilon), epsilon)
logit_p <- log(p_safe / (1 - p_safe))

# Fit spline on logit scale 
fit_ml <- lm(logit_p ~ ns_basis)

# Get fitted values (on prob scale)
fitted_logit <- predict(fit_ml)

# Get starting coefficients on logit scale -- same as fitted_logit
starting_coefs <- coef(fit_ml)
fitted_logit2 <- as.numeric(X %*% starting_coefs) # this is the same as the predict(fit_ml), just doesn't use model obj

# Transform back from logit scale
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
), file = paste0("R:/Kelly/synergy_orderly/src/sim_trial_cohort/spline_setupMali_", selected_df, ".rds"))
