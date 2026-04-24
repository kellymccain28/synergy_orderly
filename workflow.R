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
# param_ranges <- list(
#   max_SMC_kill_rate = c(3, 25),      # days, adjust to realistic range
#   smc_lambda = c(0.1, 45),   # adjust based on your hill function
#   smc_kappa = c(0.5, 5),    # adjust based on your hill function
#   season_start_day = c(20, 120)            # days, adjust as needed
# )
# # Generate LHS samples (increase to 200-500 for better coverage)
# n_lhs <- 20
# lhs_design <- randomLHS(n_lhs, 4)
# # Scale to parameter ranges
# lhs_params <- data.frame(
#   max_SMC_kill_rate = lhs_design[,1] * diff(param_ranges$max_SMC_kill_rate) + param_ranges$max_SMC_kill_rate[1],
#   smc_lambda = lhs_design[,2] * diff(param_ranges$smc_lambda) + param_ranges$smc_lambda[1],
#   smc_kappa = lhs_design[,3] * diff(param_ranges$smc_kappa) + param_ranges$smc_kappa[1],
#   season_start_day = lhs_design[,4] * diff(param_ranges$season_start_day) + param_ranges$season_start_day[1]
# )
# saveRDS(lhs_params,file = paste0("R:/Kelly/synergy_orderly/cohort_param_sets/lhs_params", Sys.Date(), ".rds"))
# # send all of these runs to the cluster, iterating through the
# cohortruns <- task_create_bulk_expr(
#   orderly::orderly_run(name = 'sim_cohort',
#                         parameters = list(N = 3200,
#                                           trial_ts = 365*3,
#                                           burnin = 90,
#                                           max_SMC_kill_rate = max_SMC_kill_rate,
#                                           smc_lambda = smc_lambda,
#                                           smc_kappa = smc_kappa,
#                                           season_start_day = season_start_day,
#                                           sim_allow_superinfections = TRUE)),
#   lhs_params
# )
# hipercow_bundle_status('freezable_kakarikis')
# hipercow_bundle_log_value('freezable_kakarikis')[[10]]


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
nparams = 32*3
ncores = if(nparams > 32) 32 else nparams
rtss_grid_final_1000threshold_0liver_2000_vmin_lognormal_pbite015 <- task_create_expr(expr = run_grid_rtss(path = "R:/Kelly/synergy_orderly",
                                                                            N = 2000,
                                                                            threshold = 1000,
                                                                            n_param_sets = nparams),
                                       resources = hipercow_resources(cores = ncores))
task_log_show(rtss_grid_final_1000threshold_0liver_2000_vmin_lognormal_pbite015)

# run generic cohort 
# hipercow_environment_create(name = 'generic',
#                             sources = c("shared/rtss.R",
#                                         "shared/helper_functions.R",
#                                         "shared/cohort_sim_utils.R",
#                                         'src/sim_cohort_generic/sim_cohort_generic.R',
#                                         'src/sim_trial_cohort/sim_trial_cohort.R'
# ))
nparams = 32*2
ncores = if(nparams > (32-4)) 32 else nparams + 4
task_create_expr(sim_cohort_generic(trial_ts = 365*3,
                                    treatment_prob = 0.9, # default is 0.9 (which gives children prophylaxis)
                                    season_start_day = 147, # 137 is to start on August 15 (days since April 1) -- this was good for highly seasonal;
                                    # 115 days is July 25 which is close to trial dates; 122 is Aug 1
                                    vax_day = 30, #75;# for perennial, trying an earlier start day of vaccination so they overlap less
                                    threshold = 5000, # default is 5000 parasites per microL
                                    country_to_run = 'generic',
                                    season = 'seasonal',
                                    N = 600,
                                    n_param_sets = nparams,
                                    get_parasit = TRUE,
                                    notes = 'season from 147 and vax from 30 (scen 5); 365*3; threshold at 5000, log normal weighting, seasonal generic, fitted rtss to pbite = 0.004 (1.277, 6.203, 0.0028) and smc pars (2.37, 18.5, 0.337)'),
                 environment = 'generic',
                 resources = hipercow_resources(cores = ncores))
source('R:/Kelly/synergy_orderly/src/sim_cohort_generic/extract_sim_notes.R')

# send jobs to run generic simulation to cluster by parameter dataset: 
# pars <- readxl::read_xlsx("R:/Kelly/synergy_orderly/figures/parameters_tbl.xlsx",
#                  sheet = 'Sheet4')

# get a range of timings for vaccine and smc 
earliest_vax_3rd <- as.Date('2017-05-01')
latest_vax_3rd <- as.Date('2017-07-01')
earliest_smc_start <- as.Date('2017-07-01')
latest_smc_start <- as.Date('2017-09-01')

vax_day <- seq(earliest_vax_3rd, latest_vax_3rd, by = '2 weeks')# '2 weeks')
vax_day <- vax_day - as.Date('2017-04-01')
season_start_day <- seq(earliest_smc_start, latest_smc_start, by = '2 weeks')#'2 weeks')
season_start_day <- season_start_day - as.Date('2017-04-01')

combos <- crossing(vax_day, season_start_day)
combos$season = 'seasonal'
combos$get_parasit = FALSE
combos$notes = paste(combos$vax_day, combos$season_start_day, combos$season, sep = '_')
combos <- combos %>%
  mutate(notes = paste0('testing with gen adaptive only; ', notes))#paste0('testing with no rtss decay; ', notes))

nparams = 32
ncores = if(nparams > (32-4)) 32 else nparams + 4

task_ids2 <- pmap(combos, function(vax_day, season_start_day, get_parasit, season, notes) {
  if(get_parasit == TRUE){
    trial_ts_ = 365
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
      notes = notes
    ),
    environment = 'generic',
    resources = hipercow_resources(cores = ncores)
  )
})

# Monitor progress
task_status(unlist(task_ids))
task_log_show(unlist(task_ids)[1])
task_status(unlist(task_ids2))
task_log_show(unlist(task_ids2)[1])

# run trial cohort simulation 
# hipercow_environment_create(name = 'trial_simulations',
#                             sources = c("shared/rtss.R",
#                                         "shared/helper_functions.R",
#                                         "shared/cohort_sim_utils.R",
#                                         'src/sim_trial_cohort/sim_trial_cohort.R',
#                                         'src/sim_trial_cohort/compare_incidence.R'
#                             ))
# hipercow_provision(method = 'pkgdepends', refs = c('cyphr', 'mrc-ide/hipercow@mrc-6733'),
#                    environment = 'trial_simulations')
nparams = 32
ncores = if(nparams > 32) 32 else nparams
country = 'Mali'#'BF'#
bestreps_or_syntest = 'syntest'
task_create_expr(sim_trial_cohort(trial_ts = 365*3, 
                                  treatment_prob = 1, # default is 1 (which gives children prophylaxis)
                                  threshold = 5000, # default is 5000 parasites per microL
                                  country_to_run = country, # should be BF or Mali
                                  n_param_sets = nparams,
                                  get_parasit = FALSE,
                                  path = "R:/Kelly/synergy_orderly/",
                                  bestreps_or_syntest = bestreps_or_syntest,
                                  notes = paste0(country, '; ', bestreps_or_syntest, '; best spline fit from 03-25_Mali_2 (rep 182) or 03-25_BF_2 (rep 153)')),
                 environment = 'trial_simulations',
                 resources = hipercow_resources(cores = ncores))
# task_log_show('0c4f5a7e7f234a50289202c7200254bb') # Mali 2-arm (BEST 182)
# task_log_show('1324bce5c96cb7457c883ab715d95104') # BF 2-arm (BEST 153)
# task_log_show('efff9ae86637ca336d16382c82a6f9d2') # Mali 2-arm (SECOND BEST 138)
# task_log_show('af48c9bfc36a530ad02b8caa0a630903') # BF 2-arm (SECOND BEST 199)
# task_log_show('eb1d621f210964ed41b1e5cf6792b87b') # Mali 2-arm (BEST 182) syntest
# task_log_show('c3ae9f2ecd024450b2369fe69c0ad704') # BF 2-arm (BEST 153) syntest
# task_log_show('4591980ec00b4b26729176d5b7a848f7') # Mali 2-arm (SECOND BEST 138) syntest
# task_log_show('b025459ea2d47d835b921b5b92e618fd') # BF 2-arm (SECOND BEST 199) syntest


# Fitting the spline for chapter 6 
# hipercow_environment_create(name = 'trial_fitting',
#                             sources = c("shared/rtss.R",
#                                         "shared/helper_functions.R",
#                                         "shared/cohort_sim_utils.R",
#                                         'src/sim_trial_cohort/sim_trial_cohort.R',
#                                         'src/sim_trial_cohort/compare_incidence.R',
#                                         'src/sim_trial_cohort/optimise_sim_trial_cohort.R',
#                                         'src/sim_trial_cohort/optimisation_helpers.R'
#                             ))
# hipercow_provision(method = 'pkgdepends', refs = c('cyphr', 'mrc-ide/hipercow@mrc-6733'),
#                    environment = 'trial_fitting')
nparams = 1
ncores = 1
country = 'Mali'
task_create_expr(optimise_sim_trial_cohort(trial_ts = 365*3, 
                                           country_to_run = country, # should be BF or Mali
                                           threshold = 5000,
                                           n_param_sets = nparams,
                                           arms_to_fit = c('rtss','smc','both'),
                                           notes = paste0(country, ', optimisation to 3 arms; 200 max iterations, trying spline w 13 knots; outputting rmse over time and/or by arm; lower tol (1e-4)')),
                 environment = 'trial_simulations',
                 resources = hipercow_resources(cores = ncores))
# task_log_show('b122962afd4b05a694e47e39a36ad566') # BF 150
# task_log_show('ae6ed35b42db01ca4058a706255ad1f5') # Mali 150
# task_log_show('90f43a6659d5a183511421133f18acfa') # BF 200
# task_log_show('98ca64f0d1d87f06e0ffd38c6017024a') # Mali 200
task_log_show('bdb333c9db09a51471f4ba395e45af40') # BF 200, 14 knots 
task_log_show('5a440e11cc95a539e8d23cba1d75d42a') # Mali 200, 13 knots 
task_log_show('a56e2e66fba175fa4dd4f7f673ce25fc') # BF 200, 14 knots, rtss and smc only not both 
task_log_show('39d97b44589e1ff3a2a645a25b0a0fcf') # Mali 200, 13 knots, rtss and smc only not both 


# # compare model_trial 
# compare_task <- task_create_expr(orderly::orderly_run(name = 'compare_model_trial'))
# task_log_show(compare_task)


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

