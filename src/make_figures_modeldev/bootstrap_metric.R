# Function to bootstrap one time point
bootstrap_ratio_ci <- function(ratios, n_boot = 1000) {
  boot_medians <- replicate(n_boot, {
    median(sample(ratios, replace = TRUE), na.rm = TRUE)
  })
  
  c(
    quantile(boot_medians, 0.025, na.rm = TRUE),
    quantile(boot_medians, 0.975, na.rm = TRUE)
  )
}

# Bootstrap function that works for any ratio/efficacy
bootstrap_metric <- function(values, transform_fn = NULL, n_boot = 10000) {
  boot_medians <- replicate(n_boot, {
    boot_sample <- sample(values, replace = TRUE)
    if (!is.null(transform_fn)) {
      boot_sample <- transform_fn(boot_sample)
    }
    median(boot_sample, na.rm = TRUE)
  })
  
  c(
    median = median(if (!is.null(transform_fn)) transform_fn(values) else values, na.rm = TRUE),
    q025 = quantile(boot_medians, 0.025, na.rm = TRUE),
    q975 = quantile(boot_medians, 0.975, na.rm = TRUE),
    q25 = quantile(boot_medians, 0.25, na.rm = TRUE),
    q75 = quantile(boot_medians, 0.75, na.rm = TRUE)
  )
}
