# Workflow to remember the order of tasks, when there is an order, and to easily run them all 

library(orderly)
library(hipercow)
library(lhs)

# Trial 
orderly_run(name = 'clean_trial_data')

orderly_run(name = 'trial_results')

# Model 
rainfall_task <- task_create_expr(orderly::orderly_run(name = 'fit_rainfall'))
task_log_show(rainfall_task)

orderly_run(name = 'run_smc_rtss',
            parameters = list(
              n_particles = 50,
              n_threads = 4,
              t_inf = 10, # time of infection relative to vaccination
              ts = ceiling(365/2), # timesteps (each one is 2 days)
              tstep = 1,
              max_SMC_kill_rate = 10,
              lambda= 39,
              kappa = 3.4,
              # SMC_decay = 0.05,
              # season_start_day = 30,
              # season_length = 150,
              # smc_interval = 30,
              inf_start = 0
            ))

# to run cohort with one set of parameters (this no longer works with process_model_output and compare as of 28 July)
param_ranges <- list(
  max_SMC_kill_rate = c(3, 25),      # days, adjust to realistic range
  smc_lambda = c(0.1, 45),   # adjust based on your hill function
  smc_kappa = c(0.5, 5),    # adjust based on your hill function
  season_start_day = c(20, 120)            # days, adjust as needed
)
# Generate LHS samples (increase to 200-500 for better coverage)
n_lhs <- 20
lhs_design <- randomLHS(n_lhs, 4)
# Scale to parameter ranges
lhs_params <- data.frame(
  max_SMC_kill_rate = lhs_design[,1] * diff(param_ranges$max_SMC_kill_rate) + param_ranges$max_SMC_kill_rate[1],
  smc_lambda = lhs_design[,2] * diff(param_ranges$smc_lambda) + param_ranges$smc_lambda[1],
  smc_kappa = lhs_design[,3] * diff(param_ranges$smc_kappa) + param_ranges$smc_kappa[1],
  season_start_day = lhs_design[,4] * diff(param_ranges$season_start_day) + param_ranges$season_start_day[1]
)
saveRDS(lhs_params,file = paste0("R:/Kelly/synergy_orderly/cohort_param_sets/lhs_params", Sys.Date(), ".rds"))
# send all of these runs to the cluster, iterating through the
cohortruns <- task_create_bulk_expr(
  orderly::orderly_run(name = 'sim_cohort',
                        parameters = list(N = 3200,
                                          trial_ts = 365*3,
                                          burnin = 90,
                                          max_SMC_kill_rate = max_SMC_kill_rate,
                                          smc_lambda = smc_lambda,
                                          smc_kappa = smc_kappa,
                                          season_start_day = season_start_day,
                                          sim_allow_superinfections = TRUE)),
  lhs_params
)
hipercow_bundle_status('freezable_kakarikis')
hipercow_bundle_log_value('freezable_kakarikis')[[10]]
# cohort_TRUE <- task_create_expr(orderly::orderly_run(name = 'sim_cohort',
#                                             parameters = list(N = 3200, # size of cohort population 
#                                                               trial_ts = 365*3,
#                                                               burnin = 90,
#                                                               max_SMC_kill_rate = 8,
#                                                               smc_lambda = 17,
#                                                               smc_kappa = 0.28,
#                                                               season_start_day = 95,
#                                                               sim_allow_superinfections = TRUE
#                                             )))
# cohort_FALSE <- task_create_expr(orderly::orderly_run(name = 'sim_cohort',
#                                              parameters = list(N = 3200, # size of cohort population 
#                                                                trial_ts = 365*3,
#                                                                burnin = 90,
#                                                                max_SMC_kill_rate = 8,
#                                                                smc_lambda = 17,
#                                                                smc_kappa = 0.28,
#                                                                season_start_day = 95,
#                                                                sim_allow_superinfections = FALSE
#                                              )))
# task_log_show(cohort_TRUE)
# task_log_show(cohort_FALSE)

# Fit SMC -- need to update the parmaeters if i want to do any lhs fitting -- atm it is just running the 'best' parameters n_param_sets times 
library(hipercow)
hipercow_provision()
hipercow_environment_create(sources = c("shared/rtss.R",
                                        "shared/helper_functions.R",
                                        "shared/cohort_sim_utils.R",
                                        "shared/likelihood.R",
                                        'src/fit_smc/fit_smc.R',
                                        'src/fit_rtss/fit_rtss.R',
                                        'src/fit_smc/run_smc_test.R',
                                        'src/sim_cohort_generic/sim_cohort_generic.R'
))
nparams = 1
ncores = if(nparams > 32) 32 else nparams
fitsmctask <- task_create_expr(expr = run_fit_smc(N = 800,
                                         n_param_sets = nparams),
                      resources = hipercow_resources(cores = ncores))
task_log_show(fitsmctask)#


n_params = 32
ncores = if(n_params > 32) 32 else n_params
test_lhs_smc_withcumul <- task_create_expr(expr = run_smc_test(N = 1000,
                                                     n_param_sets = n_params),
                                 resources = hipercow_resources(cores = ncores))
task_log_show(test_lhs_smc_withcumul)#max_SMC_kill_rate = c(1, 10),# parasites per uL per 2-day timestep lambda = c(5, 50),kappa = c(0.1, 9)
#0d15c3ef223da5a3a851ea00de6495da



# Fit RTSS
nparams = 32 # this is just the number of repetitions bc not changing any params 
ncores = if(nparams > 32) 32 else nparams
rtss50_fixedparams <- task_create_expr(expr = run_fit_rtss(path = "R:/Kelly/synergy_orderly",
                                               N = 1200,
                                               n_param_sets = nparams),
                           resources = hipercow_resources(cores = ncores))
task_log_show(rtss50_fixedparams)
# source('src/fit_smc/fit_smc.R')
# run_fit_smc(N = 100,
#             n_param_sets = 2)


# run generic cohort 
hipercow_environment_create(name = 'generic',
                            sources = c("shared/rtss.R",
                                        "shared/helper_functions.R",
                                        "shared/cohort_sim_utils.R",
                                        'src/sim_cohort_generic/sim_cohort_generic.R'
))
nparams = 20
ncores = if(nparams > 32) 32 else nparams
generic_cohort <- task_create_expr(sim_cohort_generic(trial_ts = 365*3, 
                                                      sim_allow_superinfections = FALSE, 
                                                      country_to_run = 'generic',
                                                      n_param_sets = nparams),
                                   environment = 'generic',
                                   resources = hipercow_resources(cores = ncores)) # takes 2.12 hours
task_log_show(generic_cohort)
# grid_genericF <- task_create_expr(orderly::orderly_run(name = 'sim_cohort_grid',
#                                                       parameters = list(trial_ts = 365*3,
#                                                                         sim_allow_superinfections = FALSE,
#                                                                         country_to_run = "generic",
#                                                                         n_param_sets = nparams)),
#                                 resources = hipercow_resources(cores = ncores))
# task_log_show(grid_genericF)


# run trial cohort simulation 
nparams = 50
ncores = if(nparams > 32) 32 else nparams
grid_taskbf <- task_create_expr(orderly::orderly_run(name = 'sim_cohort_grid',
                                                      parameters = list(trial_ts = 365*3,
                                                                        sim_allow_superinfections = TRUE,
                                                                        country_to_run = "BF",
                                                                        n_param_sets = nparams)),
                                resources = hipercow_resources(cores = ncores))
grid_taskmali <- task_create_expr(orderly::orderly_run(name = 'sim_cohort_grid',
                                                        parameters = list(trial_ts = 365*3,
                                                                          sim_allow_superinfections = TRUE,
                                                                          country_to_run = "Mali",
                                                                          n_param_sets = nparams)),
                                  resources = hipercow_resources(cores = ncores))
task_log_show(grid_taskbf)  #0ddd01eb14f245ba46fa7f4458593438 if 5.6 hours for 50, 1000 would take 112 hours
task_log_show(grid_taskmali)

grid_taskbf_nosuperinf <- task_create_expr(orderly::orderly_run(name = 'sim_cohort_grid',
                                                                 parameters = list(trial_ts = 365*3,
                                                                                   sim_allow_superinfections = FALSE,
                                                                                   country_to_run = "BF",
                                                                                   n_param_sets = nparams)),
                                           resources = hipercow_resources(cores = ncores))
grid_taskmali_nosuperinf <- task_create_expr(orderly::orderly_run(name = 'sim_cohort_grid',
                                                                   parameters = list(trial_ts = 365*3,
                                                                                     sim_allow_superinfections = FALSE,
                                                                                     country_to_run = "Mali",
                                                                                     n_param_sets = nparams)),
                                             resources = hipercow_resources(cores = ncores))
task_log_show(grid_taskbf_nosuperinf) #2-2.8 hours for 50 --> 50 hours for 1000 
task_log_show(grid_taskmali_nosuperinf)


# process_model_output - formatting the data to match the output from the trial 
# process_task <- task_create_expr(orderly::orderly_run(name = 'process_model_output'))
# task_log_show(process_task)
process_bfF <- task_create_expr(orderly::orderly_run(name = 'process_model_output',
                                                      parameters = list(trial_ts = 365*3,
                                                                        sim_allow_superinfections = FALSE,
                                                                        country_to_run = 'BF',
                                                                        n_param_sets = nparams)))
process_maliF <- task_create_expr(orderly::orderly_run(name = 'process_model_output',
                                                        parameters = list(trial_ts = 365*3,
                                                                          sim_allow_superinfections = FALSE,
                                                                          country_to_run = 'Mali',
                                                                          n_param_sets = nparams)))
process_bfT <- task_create_expr(orderly::orderly_run(name = 'process_model_output',
                                                      parameters = list(trial_ts = 365*3, 
                                                                        sim_allow_superinfections = TRUE,
                                                                        country_to_run = 'BF',
                                                                        n_param_sets = nparams)))
process_maliT <- task_create_expr(orderly::orderly_run(name = 'process_model_output',
                                                        parameters = list(trial_ts = 365*3, 
                                                                          sim_allow_superinfections = TRUE,
                                                                          country_to_run = 'Mali',
                                                                          n_param_sets = nparams)))
task_log_show(process_bfF)
task_log_show(process_maliF)
task_log_show(process_bfT)
task_log_show(process_maliT)


# compare model_trial 
compare_task <- task_create_expr(orderly::orderly_run(name = 'compare_model_trial'))
task_log_show(compare_task)
