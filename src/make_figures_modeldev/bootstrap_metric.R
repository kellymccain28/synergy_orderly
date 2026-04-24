# Function to bootstrap one time point
# bootstrap_ratio_ci <- function(ratios, n_boot = 5000) {
#   boot_medians <- replicate(n_boot, {
#     median(sample(ratios, replace = TRUE), na.rm = TRUE)
#   })
#   
#   c(
#     quantile(boot_medians, 0.025, na.rm = TRUE),
#     quantile(boot_medians, 0.975, na.rm = TRUE)
#   )
# }

# Bootstrap function that works for any ratio/efficacy
bootstrap_metric <- function(values, transform_fn = NULL, n_boot = 10000) {
  if (!is.null(transform_fn)) values <- transform_fn(values)
  
  boot_medians <- replicate(n_boot, {
    boot_sample <- sample(values, replace = TRUE)
    median(boot_sample, na.rm = TRUE)
  })
  
  c(
    median = median(values, na.rm = TRUE),
    q025 = quantile(boot_medians, 0.025, na.rm = TRUE),
    q975 = quantile(boot_medians, 0.975, na.rm = TRUE)
  )
}

# Boostrap IRR -- unit of resampling is individual replicate
# bootstrap_irr <- function(cases_a, pm_a, cases_b, pm_b, 
#                           n_boot = 10000, transform_fn = NULL) {
#   n <- length(cases_a)
#   boot_ratios <- replicate(n_boot, {
#     idx <- sample(n, replace = TRUE)
#     rate_a <- sum(cases_a[idx]) / sum(pm_a[idx])
#     rate_b <- sum(cases_b[idx]) / sum(pm_b[idx])
#     ratio  <- rate_a / rate_b
#     if (!is.null(transform_fn)) transform_fn(ratio) else ratio
#   })
#   c(
#     median = median(if (!is.null(transform_fn)) 
#       transform_fn(cases_a/pm_a * (sum(pm_a)/sum(cases_a))) 
#       else sum(cases_a)/sum(pm_a) / (sum(cases_b)/sum(pm_b)), 
#       na.rm = TRUE),
#     q025 = quantile(boot_ratios, 0.025, na.rm = TRUE),
#     q975 = quantile(boot_ratios, 0.975, na.rm = TRUE),
#     q25  = quantile(boot_ratios, 0.25,  na.rm = TRUE),
#     q75  = quantile(boot_ratios, 0.75,  na.rm = TRUE)
#   )
# }