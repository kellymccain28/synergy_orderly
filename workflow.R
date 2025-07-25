# Workflow to remember the order of tasks, when there is an order, and to easily run them all 

library(orderly2)

# Trial 
orderly_run(name = 'clean_trial_data')

orderly_run(name = 'trial_results')

# Model 
orderly_run(name = 'fit_rainfall')

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


orderly_run(name = 'sim_cohort',
            parameters = list(N = 2000#, # size of cohort population 
                              # max_SMC_kill_rate = 8,
                              # SMC_decay = 0.05,
                              # season_start_day = 50,#213,
                              # season_length = 30*5,
                              # smc_interval = 30
            ))

# process_model_output - formatting the data to match the output from the trial 
orderly_run(name = 'process_model_output')


# compare model_trial 
orderly_run(name = 'compare_model_trial')
