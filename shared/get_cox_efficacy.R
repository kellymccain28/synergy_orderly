# Function to get Cox PH efficacy as in Chandramohan et al. 2021
# by year and overall 

# df -- formatted dataset of person time
get_cox_efficacy <- function(df, 
                             ref,
                             model = FALSE){
  
  tidy_results <- vector("list", 4)

  # Make different columns for various ref = 
  df <- df %>%
    mutate(arm_smcref = factor(arm, levels = c('smc','rtss','both','none')),
           arm_rtssref = factor(arm, levels = c('rtss','smc','both','none')),
           arm_noneref = factor(arm, levels = c('none','rtss','smc','both')),
           end_time = ifelse(start_time==end_time, end_time + 0.00001, end_time))
           
  for (i in 1:4){
    
    if(i < 4){
      d <- df %>% 
        dplyr::filter(year == i)
    } else d <- df
    
    # Create formula dynamically using the ref argument
    if(model){
      formula_str <- paste("Surv(start_time, end_time, event) ~ factor(", ref, ")")
      cox_formula <- as.formula(formula_str)
      
      clustervar = d$rid#child_id
    } else {
      if(length(unique(df$sim_id)) >= 1){
        formula_str <- paste("Surv(start_time, end_time, event) ~ factor(", ref, ") + factor(country)")
      } else {
        formula_str <- paste("Surv(start_time, end_time, event) ~ factor(", ref, ")")
      }
      cox_formula <- as.formula(formula_str)
      
      clustervar = d$rid
    }
    
    coxmodel <- coxph(cox_formula,# + factor(country), 
                                   data = d,
                                   cluster = clustervar,
                                   ties = "efron")
    
    results <- tidy(coxmodel,
                    exponentiate = TRUE,
                    conf.int = TRUE) %>%
      # Calculate vaccine efficacy and its confidence intervals
      mutate(
        term = case_when(
          term == 'factor(arm_rtssref)smc' ~ "SMC vs. RTSS",
          term == 'factor(arm_rtssref)both' ~ "Both vs. RTSS",
          term == 'factor(arm_rtssref)vaxsmc' ~ "Both vs. RTSS",
          term == 'factor(arm_rtssref)none' ~ "None vs. RTSS",
          term == 'factor(country)Mali' ~ 'Mali vs. BF',
          term == 'factor(arm_smcref)rtss' ~ "RTSS vs. SMC",
          term == 'factor(arm_smcref)both' ~ "Both vs. SMC",
          term == 'factor(arm_smcref)vax' ~ "RTSS vs. SMC",
          term == 'factor(arm_smcref)vaxsmc' ~ "Both vs. SMC",
          term == 'factor(arm_smcref)none' ~ "None vs. SMC",
          term == 'factor(arm_noneref)vax' ~ "RTSS vs None",
          term == 'factor(arm_noneref)rtss' ~ "RTSS vs None",
          term == 'factor(arm_noneref)smc' ~ "SMC vs None",
          term == 'factor(arm_noneref)vaxsmc' ~ "Both vs None",
          term == 'factor(arm_noneref)both' ~ "Both vs None",
          TRUE ~ term
        ),
        VE = 1 - estimate,              # VE = 1 - HR
        VE_lower = 1 - conf.high,       # Lower CI for VE = 1 - Upper CI for HR
        VE_upper = 1 - conf.low ,        # Upper CI for VE = 1 - Lower CI for HR
        year = ifelse(i == 4, 'overall',as.character(i))
      ) %>% 
      rename(
        HR = estimate, 
        HR_lower = conf.low, 
        HR_upper = conf.high,
      ) %>%
      mutate(
        n_events = coxmodel$nevent,
        n_obs = coxmodel$n
      )
    
    tidy_results[[i]] <- results
  }
  
  tidy_results <- bind_rows(tidy_results)
  
  return(tidy_results)
}


# function to split by both scenario AND repetition
# get_cox_efficacy_by_repetition <- function(df, ref, model = FALSE) {
#   
#   # Get unique repetitions (sim_id in your data)
#   reps <- unique(df$sim_id)
#   
#   # Make different columns for various ref = 
#   df <- df %>%
#     mutate(arm_smcref = factor(arm, levels = c('smc','rtss','both','none')),
#            arm_rtssref = factor(arm, levels = c('rtss','smc','both','none')),
#            arm_noneref = factor(arm, levels = c('none','rtss','smc','both')),
#            end_time = ifelse(start_time==end_time, end_time + 0.00001, end_time))
#   
#   # Store results for all repetitions
#   all_results <- list()
#   
#   for (r in reps) {
#     
#     # Filter to one repetition
#     d_rep <- df %>% filter(sim_id == r)
#     
#     # Now process each year (1, 2, 3, overall)
#     yearly_results <- vector("list", 4)
#     
#     for (i in 1:4) {
#       
#       if(i < 4) {
#         d <- d_rep %>% filter(year == i)
#       } else {
#         d <- d_rep  # overall
#       }
#       
#       # Skip if no data
#       if(nrow(d) == 0) next
#       
#       # Create formula dynamically using the ref argument
#       if(model) {
#         # For model = TRUE, we have multiple sims? No, we're already in one sim_id
#         # But keep the pattern from original
#         if(length(unique(d_rep$sim_id)) == 1) {
#           formula_str <- paste("Surv(start_time, end_time, event) ~ factor(", ref, ")")
#         } else {
#           formula_str <- paste("Surv(start_time, end_time, event) ~ factor(", ref, ") + strata(sim_id)")
#         }
#         cox_formula <- as.formula(formula_str)
#         clustervar <- d$rid  # using rid as in original
#       } else {
#         if(length(unique(d_rep$sim_id)) >= 1) {
#           formula_str <- paste("Surv(start_time, end_time, event) ~ factor(", ref, ") + factor(country)")
#         } else {
#           formula_str <- paste("Surv(start_time, end_time, event) ~ factor(", ref, ")")
#         }
#         cox_formula <- as.formula(formula_str)
#         clustervar <- d$rid
#       }
#       
#       # Fit Cox model
#       coxmodel <- tryCatch({
#         coxph(cox_formula, 
#               data = d,
#               cluster = clustervar,
#               ties = "efron")
#       }, error = function(e) {
#         warning(paste("Cox failed for sim_id", r, "year", i, ":", e$message))
#         return(NULL)
#       })
#       
#       if(is.null(coxmodel)) next
#       
#       # Extract results (matching original output format)
#       results <- tidy(coxmodel,
#                       exponentiate = TRUE,
#                       conf.int = TRUE) %>%
#         mutate(
#           term = case_when(
#             term == 'factor(arm_rtssref)smc' ~ "SMC vs. RTSS",
#             term == 'factor(arm_rtssref)both' ~ "Both vs. RTSS",
#             term == 'factor(arm_rtssref)vaxsmc' ~ "Both vs. RTSS",
#             term == 'factor(arm_rtssref)none' ~ "None vs. RTSS",
#             term == 'factor(country)Mali' ~ 'Mali vs. BF',
#             term == 'factor(arm_smcref)rtss' ~ "RTSS vs. SMC",
#             term == 'factor(arm_smcref)both' ~ "Both vs. SMC",
#             term == 'factor(arm_smcref)vax' ~ "RTSS vs. SMC",
#             term == 'factor(arm_smcref)vaxsmc' ~ "Both vs. SMC",
#             term == 'factor(arm_smcref)none' ~ "None vs. SMC",
#             term == 'factor(arm_noneref)vax' ~ "RTSS vs None",
#             term == 'factor(arm_noneref)rtss' ~ "RTSS vs None",
#             term == 'factor(arm_noneref)smc' ~ "SMC vs None",
#             term == 'factor(arm_noneref)vaxsmc' ~ "Both vs None",
#             term == 'factor(arm_noneref)both' ~ "Both vs None",
#             TRUE ~ term
#           ),
#           VE = 1 - estimate,
#           VE_lower = 1 - conf.high,
#           VE_upper = 1 - conf.low,
#           year = ifelse(i == 4, 'overall', as.character(i))
#         ) %>%
#         rename(
#           HR = estimate,
#           HR_lower = conf.low,
#           HR_upper = conf.high,
#         ) %>%
#         mutate(
#           n_events = coxmodel$nevent,
#           n_obs = coxmodel$n,
#           sim_id = r  # Add the simulation ID
#         )
#       
#       yearly_results[[i]] <- results
#     }
#     
#     all_results[[r]] <- bind_rows(yearly_results)
#   }
#   
#   # Combine all repetitions
#   final_results <- bind_rows(all_results)
#   
#   return(final_results)
# }

# to get the monthly efficacy, i think i'll need to reformat the data to be grouped by month, to have person time 
# in each month over the whole study period -- or since we are assuming no loss of follow up here, could just use 
# n as person time and have the number of cases per month / personmonths (just n)?

