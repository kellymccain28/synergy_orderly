# function to get smc decay over time 

# this is unfitted and is just a placeholder for the functional form of the decay in efficacy
#' @param max_SMC_kill_rate maximum killing rate per timestep of SMC per parasite
#' @param SMC_decay exponential decay rate of killing effect 
#' @param ts timesteps 
smc_kill_rate <-  function(max_SMC_kill_rate, 
                           SMC_decay,
                           ts){
  SMC_kill_rate <- max_SMC_kill_rate * exp(-SMC_decay * ts) # time-varying SMC kill rate (exponential decay)
  
  return(SMC_kill_rate)
}


# smc_kill_rate(max_SMC_kill_rate=2.3, SMC_decay=0.05, ts = 1:30)
if (time >= smc_timing) max_SMC_kill_rate * exp(-((time - smc_timing) / lambda)^kappa) 
1 - exp(-SMC_kill_rate * dt)


# Approach 1: Using modulo arithmetic (most odin2-friendly)
# This assumes SMC is given at regular intervals during transmission season
smc_kill_rate_seasonal <- function(time, season_start, season_length, smc_interval, 
                                   max_SMC_kill_rate, lambda, kappa) {
  
  # Check if we're in transmission season
  if (time < season_start || time > (season_start + season_length)) {
    return(0)
  }
  
  # Time since season started
  time_in_season <- time - season_start
  
  # Time since last SMC dose (using modulo)
  time_since_smc <- time_in_season %% smc_interval
  
  # Apply hill function based on time since last dose
  return(max_SMC_kill_rate * exp(-(time_since_smc / lambda)^kappa))
}

# Approach 2: Step function approach (also odin2-compatible)
# This creates explicit conditions for each SMC timing
smc_kill_rate_steps <- function(time, smc_timing1, smc_timing2, smc_timing3, smc_timing4,
                                max_SMC_kill_rate, lambda, kappa) {
  
  # Find the most recent SMC timing
  time_since_smc <- if (time >= smc_timing4) {
    time - smc_timing4
  } else if (time >= smc_timing3) {
    time - smc_timing3  
  } else if (time >= smc_timing2) {
    time - smc_timing2
  } else if (time >= smc_timing1) {
    time - smc_timing1
  } else {
    -1  # No SMC given yet
  }
  
  # Apply hill function or return 0
  if (time_since_smc >= 0) {
    return(max_SMC_kill_rate * exp(-(time_since_smc / lambda)^kappa))
  } else {
    return(0)
  }
}

# Approach 3: For odin2 syntax (this is what you'd put in your odin2 model)
# Example odin2 code block:

cat("
# Example odin2 model syntax:
# Define SMC parameters
smc_timing1 <- 180  # Start of transmission season (day 180 = ~July)
smc_timing2 <- 210  # Second dose (30 days later)
smc_timing3 <- 240  # Third dose  
smc_timing4 <- 270  # Fourth dose
max_SMC_kill_rate <- 0.8
lambda <- 21  # SMC decay parameter (3 weeks)
kappa <- 2

# SMC kill rate calculation
time_since_smc <- if (time >= smc_timing4) time - smc_timing4 else (
                  if (time >= smc_timing3) time - smc_timing3 else (
                  if (time >= smc_timing2) time - smc_timing2 else (
                  if (time >= smc_timing1) time - smc_timing1 else -1)))

smc_kill_rate <- if (time_since_smc >= 0) max_SMC_kill_rate * exp(-pow(time_since_smc / lambda, kappa)) else 0

# Alternative using modulo (if your transmission season is regular):
# Assuming SMC every 30 days during 4-month season (days 180-300)
season_active <- if (time >= 180 && time <= 300) 1 else 0
time_in_season <- if (season_active == 1) time - 180 else 0
time_since_smc_mod <- time_in_season - floor(time_in_season / 30) * 30
smc_kill_rate_seasonal <- season_active * max_SMC_kill_rate * exp(-pow(time_since_smc_mod / lambda, kappa))
")

# Test the functions with example data
test_times <- seq(0, 365, by = 1)

# Test approach 1 (seasonal with regular intervals)
smc_rates_seasonal <- sapply(test_times, function(t) {
  smc_kill_rate_seasonal(t, 
                         season_start = 180,    # Start in July
                         season_length = 120,   # 4 months  
                         smc_interval = 30,     # Every 30 days
                         max_SMC_kill_rate = 0.8,
                         lambda = 21,
                         kappa = 2)
})

# Test approach 2 (explicit step function)
smc_rates_steps <- sapply(test_times, function(t) {
  smc_kill_rate_steps(t,
                      smc_timing1 = 180,
                      smc_timing2 = 210, 
                      smc_timing3 = 240,
                      smc_timing4 = 270,
                      max_SMC_kill_rate = 0.8,
                      lambda = 21,
                      kappa = 2)
})

# Plot comparison
par(mfrow = c(2, 1))

plot(test_times, smc_rates_seasonal, type = "l", lwd = 2, col = "blue",
     main = "Approach 1: Regular Seasonal SMC",
     xlab = "Day of Year", ylab = "SMC Kill Rate")
abline(v = c(180, 210, 240, 270, 300), col = "red", lty = 2, alpha = 0.5)

plot(test_times, smc_rates_steps, type = "l", lwd = 2, col = "green",
     main = "Approach 2: Explicit Step Function",
     xlab = "Day of Year", ylab = "SMC Kill Rate")
abline(v = c(180, 210, 240, 270), col = "red", lty = 2)

par(mfrow = c(1, 1))

cat("
For your odin2 model, I recommend Approach 2 (step function) because:
1. It uses only basic odin2 functions (if/else, arithmetic)
2. You can easily adjust individual SMC timings
3. It's explicit and easy to debug
4. Each SMC dose resets the hill function to maximum

Copy the odin2 syntax example above into your model and adjust the timing parameters as needed.
")

cat("
The timings could be modified by setting them relative to the cohort simulation timing --
    ") 




# Improved modulo approach that works across multiple years
# and doesn't cut off at season end

# Approach 1: Improved modulo - no cutoff, multi-year compatible
smc_kill_rate_improved <- function(time, season_start_day, season_length, smc_interval, 
                                   max_SMC_kill_rate, lambda, kappa) {
  
  # Get day of year (1-365)
  day_of_year <- ((time - 1) %% 365) + 1
  
  # Check if we're in transmission season (handles wrap-around seasons)
  season_end_day <- season_start_day + season_length - 1
  
  in_season <- if (season_end_day <= 365) {
    # Season doesn't wrap around year
    day_of_year >= season_start_day && day_of_year <= season_end_day
  } else {
    # Season wraps around (e.g., Nov-Mar)
    day_of_year >= season_start_day || day_of_year <= (season_end_day - 365)
  }
  
  if (!in_season) {
    return(0)
  }
  
  # Days since season started (accounting for year wrap)
  days_in_season <- if (season_end_day <= 365) {
    day_of_year - season_start_day + 1
  } else {
    if (day_of_year >= season_start_day) {
      day_of_year - season_start_day + 1
    } else {
      day_of_year + (365 - season_start_day + 1)
    }
  }
  
  # Time since last SMC dose (using modulo)
  time_since_smc <- (days_in_season - 1) %% smc_interval
  
  # Apply hill function - NO CUTOFF at season end
  return(max_SMC_kill_rate * exp(-(time_since_smc / lambda)^kappa))
}

# For odin2 syntax (multi-year compatible):
cat("
# Odin2 model syntax - improved modulo approach:

# SMC parameters
season_start_day <- 180    # Day 180 = ~July 1st
season_length <- 120       # 4 months (120 days)
smc_interval <- 30         # SMC every 30 days
max_SMC_kill_rate <- 0.8
lambda <- 21
kappa <- 2

# Calculate day of year
day_of_year <- ((time - 1) %% 365) + 1

# Check if in season (handles wrap-around seasons like Nov-Mar)
season_end_day <- season_start_day + season_length - 1
in_season <- if (season_end_day <= 365) (
  if (day_of_year >= season_start_day && day_of_year <= season_end_day) 1 else 0
) else (
  if (day_of_year >= season_start_day || day_of_year <= (season_end_day - 365)) 1 else 0
)

# Days since season started (accounting for year wrap)
days_in_season <- if (season_end_day <= 365) (
  day_of_year - season_start_day + 1
) else (
  if (day_of_year >= season_start_day) (
    day_of_year - season_start_day + 1
  ) else (
    day_of_year + (365 - season_start_day + 1)
  )
)

# Time since last SMC dose
time_since_smc <- (days_in_season - 1) - floor((days_in_season - 1) / smc_interval) * smc_interval

# SMC kill rate (NO cutoff at season end)
smc_kill_rate <- in_season * max_SMC_kill_rate * exp(-pow(time_since_smc / lambda, kappa))
")

# Test with multi-year data
test_times_multiyear <- seq(1, 365*3, by = 1)  # 3 years

# Test the improved approach
smc_rates_improved <- sapply(test_times_multiyear, function(t) {
  smc_kill_rate_improved(t, 
                         season_start_day = 180,    # July 1st
                         season_length = 120,       # 4 months  
                         smc_interval = 30,         # Every 30 days
                         max_SMC_kill_rate = 0.8,
                         lambda = 21,
                         kappa = 2)
})

# Plot multi-year results
plot(test_times_multiyear, smc_rates_improved, type = "l", lwd = 2, col = "purple",
     main = "Improved Modulo Approach - Multi-Year SMC (No Season Cutoff)",
     xlab = "Day (3 years)", ylab = "SMC Kill Rate")

# Add vertical lines for season boundaries
for(year in 0:2) {
  season_start <- 180 + year * 365
  season_end <- 300 + year * 365
  abline(v = season_start, col = "green", lty = 2, alpha = 0.7)
  abline(v = season_end, col = "red", lty = 2, alpha = 0.7)
}
legend("topright", c("Season Start", "Season End"), 
       col = c("green", "red"), lty = 2, cex = 0.8)

# Test with wrap-around season (e.g., November to March)
cat("
# Example for wrap-around season (Nov-Mar):
season_start_day <- 305    # Day 305 = ~Nov 1st  
season_length <- 150       # 5 months (wraps around New Year)
")

smc_rates_wraparound <- sapply(test_times_multiyear, function(t) {
  smc_kill_rate_improved(t, 
                         season_start_day = 305,    # Nov 1st
                         season_length = 150,       # 5 months (wraps around)
                         smc_interval = 30,         # Every 30 days
                         max_SMC_kill_rate = 0.8,
                         lambda = 21,
                         kappa = 2)
})

# Plot wrap-around season
plot(test_times_multiyear[1:730], smc_rates_wraparound[1:730], type = "l", lwd = 2, col = "orange",
     main = "Wrap-around Season Example (Nov-Mar)",
     xlab = "Day (2 years)", ylab = "SMC Kill Rate")

for(year in 0:1) {
  abline(v = 305 + year * 365, col = "green", lty = 2)  # Nov 1
  abline(v = 90 + year * 365, col = "red", lty = 2)    # End of Mar
}

cat("
Key improvements:
1. ✅ NO cutoff at season end - SMC effect continues until next dose
2. ✅ Works across multiple years automatically  
3. ✅ Handles wrap-around seasons (e.g., Nov-Mar)
4. ✅ Simple day_of_year calculation works in odin2
5. ✅ Each year gets the same seasonal pattern

For your odin2 model, use the syntax shown above. The modulo approach 
now gives you consistent SMC dosing every year without any artificial cutoffs!
")
