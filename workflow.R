# Workflow to remember the order of tasks, when there is an order, and to easily run them all 

library(orderly2)
library(hipercow)

# Trial 
orderly_run(name = 'clean_trial_data')

orderly_run(name = 'trial_results')

# Model 
rainfall_task <- task_create_expr(orderly2::orderly_run(name = 'fit_rainfall'))
task_log_show(rainfall_task)

orderly_run(name = 'run_smc_rtss',
            parameters = list(
              n_particles = 50,
              n_threads = 4,
              t_inf = 10, # time of infection relative to vaccination
              ts = 365, # timesteps (each one is 2 days)
              tstep = 1,
              max_SMC_kill_rate = 10,
              SMC_decay = 0.05,
              season_start_day = 30,
              season_length = 150,
              smc_interval = 30,
              inf_start = 0
            ))

# to run cohort with one set of parameters (this no longer works with process_model_output and compare as of 28 July)
orderly_run(name = 'sim_cohort',
            parameters = list(N = 2000, # size of cohort population 
                              trial_ts = 365,
                              burnin = 50,
                              max_SMC_kill_rate = 8,
                              smc_lambda = 17,
                              smc_kappa = 0.28,
                              season_start_day = 95,
                              sim_allow_superinfections = TRUE
                              # SMC_decay = 0.05,
                              # season_start_day = 50,#213,
                              # season_length = 30*5,
                              # smc_interval = 30
            ))


grid_taskbf <- task_create_expr(orderly2::orderly_run(name = 'sim_cohort_grid',
                                                      parameters = list(trial_ts = 365*3,
                                                                        sim_allow_superinfections = TRUE,
                                                                        country_to_run = "BF",
                                                                        n_param_sets = 30)))
grid_taskmali <- task_create_expr(orderly2::orderly_run(name = 'sim_cohort_grid',
                                                        parameters = list(trial_ts = 365*3,
                                                                          sim_allow_superinfections = TRUE,
                                                                          country_to_run = "Mali",
                                                                          n_param_sets = 30)))
task_log_show(grid_taskbf) 
task_log_show(grid_taskmali)

grid_taskbf_nosuperinf <- task_create_expr(orderly2::orderly_run(name = 'sim_cohort_grid',
                                                                 parameters = list(trial_ts = 365*3,
                                                                                   sim_allow_superinfections = FALSE,
                                                                                   country_to_run = "BF",
                                                                                   n_param_sets = 30)))
grid_taskmali_nosuperinf <- task_create_expr(orderly2::orderly_run(name = 'sim_cohort_grid',
                                                                   parameters = list(trial_ts = 365*3,
                                                                                     sim_allow_superinfections = FALSE,
                                                                                     country_to_run = "Mali",
                                                                                     n_param_sets = 30)))
task_log_show(grid_taskbf_nosuperinf) 
task_log_show(grid_taskmali_nosuperinf)

# process_model_output - formatting the data to match the output from the trial 
# process_task <- task_create_expr(orderly2::orderly_run(name = 'process_model_output'))
# task_log_show(process_task)
process_bfF <- task_create_expr(orderly2::orderly_run(name = 'process_model_output',
                                                      parameters = list(trial_ts = 365*3,
                                                                        sim_allow_superinfections = FALSE,
                                                                        country_to_run = 'BF',
                                                                        n_param_sets = 30)))
process_maliF <- task_create_expr(orderly2::orderly_run(name = 'process_model_output',
                                                        parameters = list(trial_ts = 365*3,
                                                                          sim_allow_superinfections = FALSE,
                                                                          country_to_run = 'Mali',
                                                                          n_param_sets = 30)))
process_bfT <- task_create_expr(orderly2::orderly_run(name = 'process_model_output',
                                                      parameters = list(trial_ts = 365*3, 
                                                                        sim_allow_superinfections = TRUE,
                                                                        country_to_run = 'BF',
                                                                        n_param_sets = 30)))
process_maliT <- task_create_expr(orderly2::orderly_run(name = 'process_model_output',
                                                        parameters = list(trial_ts = 365*3, 
                                                                          sim_allow_superinfections = TRUE,
                                                                          country_to_run = 'Mali',
                                                                          n_param_sets = 30)))
task_log_show(process_bfF)
task_log_show(process_bfT)
task_log_show(process_maliF)
task_log_show(process_maliT)


# compare model_trial 
compare_task <- task_create_expr(orderly2::orderly_run(name = 'compare_model_trial'))
task_log_show(compare_task)
