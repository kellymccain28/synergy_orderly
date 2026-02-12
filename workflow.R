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
                                        'src/sim_cohort_generic/sim_cohort_generic.R',
                                        'src/fit_smc/calculate_efficacy_likelihood.R',
                                        'src/fit_rtss/run_grid_rtss.R'
))
# nparams = 1
# ncores = if(nparams > 32) 32 else nparams
# fitsmctask <- task_create_expr(expr = run_fit_smc(N = 2400,
#                                          n_param_sets = nparams),
#                       resources = hipercow_resources(cores = ncores))
# task_log_show(fitsmctask)#



n_params = 32*3
ncores = if(n_params > 32) 32 else n_params
test_lhs_smc_7day <- task_create_expr(expr = run_smc_test(N = 2000,
                                                     n_param_sets = n_params,
                                                     threshold = 3000), # threshold of 3000 matches the one used for clincial malaria in Zongo study 
                                 resources = hipercow_resources(cores = ncores))
task_log_show(test_lhs_smc_7day)




# Fit RTSS
# nparams = 3 # this is just the number of repetitions bc not changing any params 
# ncores = if(nparams > 32) 32 else nparams
# rtss20_fixedparams_months <- task_create_expr(expr = run_fit_rtss(path = "R:/Kelly/synergy_orderly",
#                                                N = 5000,
#                                                n_param_sets = nparams),
#                            resources = hipercow_resources(cores = ncores))
# task_log_show(rtss20_fixedparams_months)

# Grid RTSS 
nparams = 32
ncores = if(nparams > 32) 32 else nparams
rtss_grid_final_1000threshold_0liver_2000_vmin_lognormal <- task_create_expr(expr = run_grid_rtss(path = "R:/Kelly/synergy_orderly",
                                                                            N = 2000,
                                                                            threshold = 1000,
                                                                            n_param_sets = nparams),
                                       resources = hipercow_resources(cores = ncores))
task_log_show(rtss_grid_final_1000threshold_0liver_2000_vmin_lognormal)

# run generic cohort 
hipercow_environment_create(name = 'generic',
                            sources = c("shared/rtss.R",
                                        "shared/helper_functions.R",
                                        "shared/cohort_sim_utils.R",
                                        'src/sim_cohort_generic/sim_cohort_generic.R',
                                        'src/sim_trial_cohort/sim_trial_cohort.R'
))
nparams = 32*3
ncores = if(nparams > (32-4)) 32 else nparams + 4
task_create_expr(sim_cohort_generic(trial_ts = 365*3, 
                                    treatment_prob = 0.9, # default is 0.9 (which gives children prophylaxis)
                                    season_start_day = 122, # 137 is to start on August 15 (days since April 1) -- this was good for highly seasonal; 
                                    # 115 days is July 25 which is close to trial dates; 122 is Aug 1
                                    vax_day = 68, #75;# for perennial, trying an earlier start day of vaccination so they overlap less 
                                    threshold = 5000, # default is 5000 parasites per microL
                                    country_to_run = 'generic',
                                    season = 'seasonal',         
                                    N = 3000,
                                    n_param_sets = nparams,
                                    get_parasit = FALSE,
                                    notes = 'season from 122 and vax from 68; 365*3; threshold at 5000, log normal weighting, seasonal generic, fitted rtss (1.77, 2.63, 0.000513) and smc pars (2.37, 18.5, 0.337)'),
                 environment = 'generic',
                 resources = hipercow_resources(cores = ncores)) 
# 'e5eb6c04b305ab4e9695e353fc4392e1' - 
task_log_show('e5eb6c04b305ab4e9695e353fc4392e1') 
source('R:/Kelly/synergy_orderly/src/sim_cohort_generic/extract_sim_notes.R')

# send jobs to run generic simulation to cluster by parameter dataset: 
pars <- readxl::read_xlsx("R:/Kelly/synergy_orderly/figures/parameters_tbl.xlsx",
                 sheet = 'Sheet4')
nparams = 32*3
ncores = if(nparams > (32-4)) 32 else nparams + 4
task_ids <- pmap(pars, function(vax_day, season_start_day, get_parasit, season, notes) {
  if(get_parasit == TRUE){
    trial_ts_ = 50
  } else {trial_ts_ = 365 * 3}
  
  task_create_expr(
    sim_cohort_generic(
      trial_ts = trial_ts_,
      treatment_prob = 0.9,
      season_start_day = !!season_start_day,
      vax_day = !!vax_day,
      threshold = 5000,
      country_to_run = 'generic',
      season = !!season,
      N = 3000,
      n_param_sets = nparams,
      get_parasit = !!get_parasit,
      notes = paste0('season from ', !!season_start_day, ' and vax from ', !!vax_day,
                     '; 365*3; threshold at 5000, log normal weighting, ',
                     !!season, ' generic, fitted rtss (1.74, 4.69, 0.00237) and
                     smc (2.37, 18.5, 0.337) pars')
    ),
    environment = 'generic',
    resources = hipercow_resources(cores = ncores)
  )
})

# Monitor progress
task_status(unlist(task_ids))
task_log_show(unlist(task_ids)[1])

# run trial cohort simulation 
hipercow_environment_create(name = 'trial_simulations',
                            sources = c("shared/rtss.R",
                                        "shared/helper_functions.R",
                                        "shared/cohort_sim_utils.R",
                                        'src/sim_trial_cohort/sim_trial_cohort.R',
                                        'src/sim_trial_cohort/compare_incidence.R'
                            ))
# hipercow_provision(method = 'pkgdepends', refs = c('cyphr', 'mrc-ide/hipercow@mrc-6733'),
#                    environment = 'trial_simulations')
nparams = 32 
ncores = if(nparams > 32) 32 else nparams
country = 'BF'
trial_sim2 <- task_create_expr(sim_trial_cohort(trial_ts = 365*3, 
                                               treatment_prob = 0.9, # default is 1 (which gives children prophylaxis)
                                               threshold = 5000, # default is 5000 parasites per microL
                                               country_to_run = country, # should be BF or Mali
                                               n_param_sets = nparams,
                                               get_parasit = FALSE,
                                               path = "R:/Kelly/synergy_orderly/",
                                               notes = 'test simulations with grid and likelihood/rsme calculations, without offset'),
                              environment = 'trial_simulations',
                              resources = hipercow_resources(cores = ncores))
task_log_show(trial_sim2)


# run trial cohort simulation -- old
# nparams = 50
# ncores = if(nparams > 32) 32 else nparams
# grid_taskbf <- task_create_expr(orderly::orderly_run(name = 'sim_cohort_grid',
#                                                       parameters = list(trial_ts = 365*3,
#                                                                         sim_allow_superinfections = TRUE,
#                                                                         country_to_run = "BF",
#                                                                         n_param_sets = nparams)),
#                                 resources = hipercow_resources(cores = ncores))
# grid_taskmali <- task_create_expr(orderly::orderly_run(name = 'sim_cohort_grid',
#                                                         parameters = list(trial_ts = 365*3,
#                                                                           sim_allow_superinfections = TRUE,
#                                                                           country_to_run = "Mali",
#                                                                           n_param_sets = nparams)),
#                                   resources = hipercow_resources(cores = ncores))
# task_log_show(grid_taskbf)  #0ddd01eb14f245ba46fa7f4458593438 if 5.6 hours for 50, 1000 would take 112 hours
# task_log_show(grid_taskmali)
# 
# grid_taskbf_nosuperinf <- task_create_expr(orderly::orderly_run(name = 'sim_cohort_grid',
#                                                                  parameters = list(trial_ts = 365*3,
#                                                                                    sim_allow_superinfections = FALSE,
#                                                                                    country_to_run = "BF",
#                                                                                    n_param_sets = nparams)),
#                                            resources = hipercow_resources(cores = ncores))
# grid_taskmali_nosuperinf <- task_create_expr(orderly::orderly_run(name = 'sim_cohort_grid',
#                                                                    parameters = list(trial_ts = 365*3,
#                                                                                      sim_allow_superinfections = FALSE,
#                                                                                      country_to_run = "Mali",
#                                                                                      n_param_sets = nparams)),
#                                              resources = hipercow_resources(cores = ncores))
# task_log_show(grid_taskbf_nosuperinf) #2-2.8 hours for 50 --> 50 hours for 1000 
# task_log_show(grid_taskmali_nosuperinf)


# process_model_output - formatting the data to match the output from the trial 
# process_task <- task_create_expr(orderly::orderly_run(name = 'process_model_output'))
# task_log_show(process_task)
# process_bfF <- task_create_expr(orderly::orderly_run(name = 'process_model_output',
#                                                       parameters = list(trial_ts = 365*3,
#                                                                         sim_allow_superinfections = FALSE,
#                                                                         country_to_run = 'BF',
#                                                                         n_param_sets = nparams)))
# process_maliF <- task_create_expr(orderly::orderly_run(name = 'process_model_output',
#                                                         parameters = list(trial_ts = 365*3,
#                                                                           sim_allow_superinfections = FALSE,
#                                                                           country_to_run = 'Mali',
#                                                                           n_param_sets = nparams)))
# process_bfT <- task_create_expr(orderly::orderly_run(name = 'process_model_output',
#                                                       parameters = list(trial_ts = 365*3, 
#                                                                         sim_allow_superinfections = TRUE,
#                                                                         country_to_run = 'BF',
#                                                                         n_param_sets = nparams)))
# process_maliT <- task_create_expr(orderly::orderly_run(name = 'process_model_output',
#                                                         parameters = list(trial_ts = 365*3, 
#                                                                           sim_allow_superinfections = TRUE,
#                                                                           country_to_run = 'Mali',
#                                                                           n_param_sets = nparams)))
# task_log_show(process_bfF)
# task_log_show(process_maliF)
# task_log_show(process_bfT)
# task_log_show(process_maliT)


# compare model_trial 
compare_task <- task_create_expr(orderly::orderly_run(name = 'compare_model_trial'))
task_log_show(compare_task)
