# Script to investigate the synergy/antagonism 
library(hipercow)
library(ggplot2)
library(zoo)
# to plot: 
# - parasitemia over time by arm 
# - immunity (each type) over time by arm 
# - test synergy/antagonism with varied immunity (e.g. fixed gen adaptive)
# - get proportion of non-zero infections that become cases in RTSS vs no int group 
# - check that SMC is more likely to clear an infection with lower starting mero

# Set up small cohort simulation ####
source("R:/Kelly/synergy_orderly/src/sim_cohort_generic/sim_cohort_generic.R", echo = TRUE)
source("R:/Kelly/synergy_orderly/src/make_figures_modeldev/summarize_IRRs.R")
source("R:/Kelly/synergy_orderly/src/make_figures_modeldev/plot_1-IRR_average.R")
source("R:/Kelly/synergy_orderly/src/make_figures_modeldev/bootstrap_metric.R")


# hipercow::task_create_expr(sim_cohort_generic(trial_ts = 365*3,
#                                               treatment_prob = 0.9, # default is 1 (which gives children prophylaxis)
#                                               season_start_day = 110, # default is 137 to start on August 15 (days since April 1)
#                                               vax_day = 75, # default is mid june (75) for 3rd dose
#                                               threshold = 5000, # default is 5000 parasites per microL
#                                               N = 800,
#                                               country_to_run = 'generic',
#                                               season = 'seasonal',
#                                               n_param_sets = 64,
#                                               path = "R:/Kelly/synergy_orderly/",
#                                               notes = 'model with time adaptive immunity tau = 30; n=800 cohort to test, 3 years; 110, 75; 64 sims; expect antagonsim',# notes are to write down the specifics of the runs
#                                               get_parasit = TRUE),
#                            environment = 'generic',
#                            resources = hipercow_resources(cores = 32))
# task_log_show('952765157eb5b7c423c9b8eeec79de81') # 64 reps of 110,75
# task_log_show('ff5d01d1aa38b5b92e0bd72d88b39f11') # 64 reps of 150,30

# Pull out data to use in the plots ####
path = 'R:/Kelly/synergy_orderly/src/sim_cohort_generic/outputs'
cohort_folder = 'sim_cohort_generic' 
# outputsfolder = 'outputs_2026-03-04_2' # basic one, smc 150, vaccine 75, n= 600; shows very slight antagonism: -7.45
# outputsfolder = 'outputs_2026-03-04_3' # basic one, same as above but with innate immunity set to 1; 19.17
# outputsfolder = 'outputs_2026-03-04_4' # basic one, same as above but with gen adaptive immunity set to 1; 36.99
outputsfolder = 'outputs_2026-03-05'   # smc 110, vaccine 75, n = 600; all normal immunity; shows antagonism  -- used this to make the figures; -45.55
# outputsfolder = 'outputs_2026-03-05_2' # smc 110, vaccine 75, n = 600; gen adapative set to 1 -- less antagonism? yes now synergy; 24.50
# outputsfolder = 'outputs_2026-03-06_2' # smc 110, vaccine 75, n = 600; gen adaptive set to 1 and also new st immunity only based on time since infection; -8.59
# outputsfolder = 'outputs_2026-03-10'   # smc 110, vaccine 75, n = 600; gen adaptive set to 1 and also new st immunity only based on time since infection, slower decay (larger tau); 15.99
outputsfolder = 'outputs_2026-03-11'    # n=1200 150, 60; overall value is 16.19
outputsfolder = 'outputs_2026-03-11_2'  # n=1200 110, 60; overall value is -10.67
outputsfolder = 'outputs_2026-03-11_3'  # n=1200 110, 75; overall value is -6.05
outputsfolder = 'outputs_2026-03-11_4'  # n=1200 150, 75; overall value is -1.14
outputsfolder = 'outputs_2026-03-11_5'  # n=1200 100, 75; overall value is -13.95
outputsfolder = 'outputs_2026-03-11_6'  # n=1200 100, 60; overall value is 1.89
outputsfolder = 'outputs_2026-03-11_7'  # n=1200 150, 50; overall value is 7.09
outputsfolder = 'outputs_2026-03-11_8'  # n=1200 150, 30; overall value is 1.85
outputsfolder = 'outputs_2026-03-11_9'  # n=1200 150, 60 part 2; overall value is -3.42
outputsfolder = 'outputs_2026-03-11_10' # n=1200 150, 60 part 3; overall value is 1.88

outputsfolder = 'outputs_2026-03-11_11' # 110, 75; 64 runs - expect antagonism
# outputsfolder = 'outputs_2026-03-11_12' # 150, 60; 64 runs - expect synergy
outputsfolder = 'outputs_2026-03-12'    # 150, 30; 64 runs - expect synergy

outputsfolder = 'outputs_2026-03-23_3' # time adaptive tau = 30; 150, 30; 64 sims; expect synergy -- -0.3284321
outputsfolder = 'outputs_2026-03-23_4' # time adaptive tau = 30; 110, 75; 64 sims; expect antagonism -- -0.3284321

make_infection_dynamics_plots <- function(outputsfolder){
  type <- if(outputsfolder %in% c('outputs_2026-03-04_2',
                                  'outputs_2026-03-12',#**
                                  'outputs_2026-03-17_2',
                                  'outputs_2026-03-20') ){
    'synergy' 
  } else if (outputsfolder %in% c('outputs_2026-03-05',
                                  'outputs_2026-03-11_11',#**
                                  'outputs_2026-03-17',
                                  'outputs_2026-03-20_2')) {
    'antagonism'
  } else if (outputsfolder == 'outputs_2026-03-10'){
      'time_immunity'
    }
  
  if(outputsfolder %in% c('outputs_2026-03-17_2')){
    parasit <- readRDS(file.path(path, outputsfolder, 'parasitemia_sml.rds'))
  } else {
    parasit <- readRDS(file.path(path, outputsfolder, 'parasitemia.rds'))
  }
  incifiles <- list.files(file.path(path, outputsfolder), pattern = "^incidence_b", full.names = TRUE)
  inci <- bind_rows(lapply(incifiles, readRDS)) %>%
    filter(date > '2017-04-01' & date < '2018-05-01') 
  # saveRDS(inci, file.path(path, outputsfolder, 'incidence.rds'))
  formattedfiles <- list.files(file.path(path, outputsfolder), pattern = '^infs', full.names = TRUE)
  formatted <- bind_rows(lapply(formattedfiles, readRDS))
  
  mycolors <- c('both' = '#E15554', 
                'none' = '#E1BC29',
                'rtss' = '#3BB273',
                'smc' = '#7768AE',
                'SMC delivery' = '#709176',
                'RTS,S delivery' = '#470024')
  
  parasit_week <- parasit %>%
    filter(parasites > 10) %>%
    mutate(week = ceiling(time/7)) %>%
    group_by(rid, arm, week) %>%
    summarise(parasites = mean(parasites))
  
  inci_agg <- inci %>%
    mutate(time_value = yearmonth, 
           time_value_num = yearmonth,
           time_unit = 'yearmonth') %>%
    split(.$sim_id) %>%
    map_dfr(~ .x %>% dplyr::select(arm, time_value, time_value_num, time_unit, 
                                   person_months, incidence_per_1000pm) %>%
              pivot_wider(
                names_from = arm,
                values_from = c(person_months, incidence_per_1000pm),
                id_cols = c(time_value, time_value_num, time_unit)
              ),
            .id = "sim_id") %>%
    mutate(rtss_none_irr = (incidence_per_1000pm_rtss / incidence_per_1000pm_none),
           smc_none_irr = (incidence_per_1000pm_smc / incidence_per_1000pm_none),
           both_none_irr = (incidence_per_1000pm_both / incidence_per_1000pm_none),
           rtss_smc_irr = (incidence_per_1000pm_rtss / incidence_per_1000pm_smc),
           both_smc_irr = (incidence_per_1000pm_both / incidence_per_1000pm_smc),
           both_rtss_irr = (incidence_per_1000pm_both / incidence_per_1000pm_rtss),
           smc_rtss_irr = (incidence_per_1000pm_smc / incidence_per_1000pm_rtss)  )%>%
    mutate(inci_averted_model = incidence_per_1000pm_none - incidence_per_1000pm_both,
           cases_averted_model = inci_averted_model * 1000, # if pop is 1000
           inci_averted_expected = incidence_per_1000pm_none - 
             (incidence_per_1000pm_rtss * incidence_per_1000pm_smc)/incidence_per_1000pm_none,
           cases_averted_expected = inci_averted_expected * 1000, # if pop is 1000
           difference_inci_averted_pred_exp = inci_averted_model - inci_averted_expected,
           difference_cases_averted_pred_exp = cases_averted_model - cases_averted_expected) 
  
  # Get the overall value of the synergy metric xxxxx this is the average monthly synergistic effect not the overall effect 
  # inci_agg %>%
  #   summarise(mean_inci_diff = mean(difference_inci_averted_pred_exp, na.rm = TRUE))
  # 30.2 for the basic one 150 and 75
  # -23.6 for the basic one 110 and 75
  
  
  inci_overall <- inci %>%
    group_by(arm, sim_id) %>%
    summarize(n_cases = sum(n_cases),
              person_months = sum(person_months)) %>%
    mutate(time_value = 'overall', 
           time_value_num = 'overall',
           time_unit = 'overall') %>%
    mutate(incidence_per_1000pm = n_cases/ person_months * 1000) %>%
    split(.$sim_id) %>%
    map_dfr(~ .x %>% dplyr::select(arm, time_value, time_value_num, time_unit, 
                                   person_months, incidence_per_1000pm) %>%
              pivot_wider(
                names_from = arm,
                values_from = c(person_months, incidence_per_1000pm),
                id_cols = c(time_value, time_value_num, time_unit)
              ),
            .id = "sim_id")%>%
    mutate(rtss_none_irr = (incidence_per_1000pm_rtss / incidence_per_1000pm_none),
           smc_none_irr = (incidence_per_1000pm_smc / incidence_per_1000pm_none),
           both_none_irr = (incidence_per_1000pm_both / incidence_per_1000pm_none),
           rtss_smc_irr = (incidence_per_1000pm_rtss / incidence_per_1000pm_smc),
           both_smc_irr = (incidence_per_1000pm_both / incidence_per_1000pm_smc),
           both_rtss_irr = (incidence_per_1000pm_both / incidence_per_1000pm_rtss),
           smc_rtss_irr = (incidence_per_1000pm_smc / incidence_per_1000pm_rtss)  )%>%
    mutate(expected_efficacy = 1 - (rtss_none_irr * smc_none_irr),
           ratio_pred_exp = (1-both_none_irr) / expected_efficacy,
           ratio_inci_rate_only = both_none_irr / (rtss_none_irr * smc_none_irr),
           inci_averted_model = incidence_per_1000pm_none - incidence_per_1000pm_both,
           cases_averted_model = inci_averted_model * 1000, # if pop is 1000
           inci_averted_expected = incidence_per_1000pm_none - 
             (incidence_per_1000pm_rtss * incidence_per_1000pm_smc)/incidence_per_1000pm_none,
           cases_averted_expected = inci_averted_expected * 1000, # if pop is 1000
           difference_inci_averted_pred_exp = inci_averted_model - inci_averted_expected,
           difference_cases_averted_pred_exp = cases_averted_model - cases_averted_expected)
  # better metric for overall synergy 
  inci_overall$difference_inci_averted_pred_exp
  mean(inci_overall$difference_inci_averted_pred_exp)
  
  # Plot parasitaemia over time by arm ####
  parmed <- parasit_week %>%
    group_by(arm, week) %>%
    summarise(mean_pb = mean(parasites),
              median_pb = median(parasites)) 
  
  # ggplot(parmed %>% filter(week < 60)) + 
  #   geom_point(aes(x = week, y = mean_pb, color = arm, linetype = 'mean')) +
  #   geom_line(aes(x = week, y = mean_pb, color = arm, linetype = 'mean')) +
  #   # geom_point(aes(x = week, y = median_pb, color = arm, linetype = 'median')) +
  #   # geom_line(aes(x = week, y = median_pb, color = arm, linetype = 'median')) +
  #   scale_y_log10() + 
  # theme_bw() + facet_wrap(~arm)
  
  # ggplot(parmed %>% filter(week < 60)) + 
  #   geom_point(aes(x = week, y = mean_pb, color = arm)) +
  #   geom_line(aes(x = week, y = mean_pb, color = arm)) +
  #   scale_y_log10() +
  #   labs(x = 'Weeks since follow-up',
  #        y = 'Mean PRBCs per \u03bcL') +
  #   theme_bw()
  
  # ggplot(parmed %>% filter(week < 60)) + 
  #   geom_point(aes(x = week, y = median_pb, color = arm)) +
  #   geom_line(aes(x = week, y = median_pb, color = arm)) +
  #   scale_y_log10() + 
  #   theme_bw() + facet_wrap(~arm)
  
  medianprbcs <- ggplot(parmed %>% filter(week < 70)) + 
    geom_point(aes(x = week, y = median_pb, color = arm), alpha = 0.9) +
    geom_line(aes(x = week, y = median_pb, color = arm), alpha=0.5) +
    scale_color_manual(values = mycolors) +
    scale_y_log10(breaks = c(1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6)) +
    scale_x_continuous(breaks = seq(0,100,10)) + 
    labs(x = 'Weeks since follow-up',
         y = 'Median PRBCs per \u03bcL') +
    theme_bw(base_size =  14)
  ggsave(paste0('R:/Kelly/synergy_orderly/src/sim_cohort_generic/plots_infection_dynamics/median_prbcs_byarm_', type,  outputsfolder,'.pdf'), medianprbcs)
  
  parasit %>% 
    group_by(arm) %>% 
    filter(parasites > 100) %>% # filter to densities above microscopy limit 10.1038/nrmicro3364 (conservative, could be 20-50) 10.1086/427243
    summarise(meanpb= mean(parasites, na.rm = TRUE),
              medianpb = median(parasites, na.rm = TRUE))
  
  # in the trial, the 'both' arm had the highest parasitemia 
  weekly <- cyphr::decrypt(readRDS('R:/Kelly/synergy_orderly/archive/clean_trial_data/20260324-104428-5bccd354/data/weekly.rds'), key) %>%
    mutate(type = 'Trial') %>%
    select(arm, pb = pf_asex_fdensity, type)
  parmed2 <- parmed %>%
    mutate(type = 'Model') %>%
    select(type, arm, pb = median_pb)
  # parmed2 <- parasit %>%
  #   ungroup() %>% filter(parasites > 10) %>%
  #   mutate(type = 'Model') %>%
  #   select(type, arm, pb = parasites)
  
  pb_model_trial <- bind_rows(weekly, parmed2) %>%
    mutate(arm = factor(arm, levels = c('both','rtss','smc','none')))
  
  pb_model_trial_plot <- ggplot(pb_model_trial) +
    geom_violin(aes(x = arm, y = pb, fill = type), 
                alpha = 0.4, position = position_dodge(0.9)) +
    geom_jitter(aes(x = arm, y = pb, color = type), 
                position = position_dodge(0.9), alpha = 0.6) +
    geom_hline(yintercept = 5000, linetype = 2) +
    scale_color_manual(values = c('Trial' = '#E15554', 
                                  'Model' = '#E1BC29')) +
    scale_fill_manual(values = c('Trial' = '#E15554', 
                                 'Model' = '#E1BC29')) +
    scale_y_log10() +
    theme_classic() +
    labs(x = NULL, color = NULL, fill = NULL,
         y = 'PRBCs per \u03bcL')
  ggsave(paste0('R:/Kelly/synergy_orderly/src/sim_cohort_generic/plots_infection_dynamics/trial_model_prbcs_byarm', type,'_',outputsfolder,'.pdf'), 
         pb_model_trial_plot, height = 6, width = 8)
  
  # ggplot(parasit %>% filter(parasites NULL# ggplot(parasit %>% filter(parasites > 100)) +
  #   geom_jitter(aes(x = arm, y = parasites)) + #, color = det, shape = country
  #   geom_violin(aes(x = arm, y= parasites), alpha = 0.4) +
  #   scale_y_log10() +
  #   theme_classic()
  # ggplot(parasit %>% filter(parasites > 0)) + 
  #   geom_jitter(aes(x = time, y = parasites, color = det)) + 
  #   geom_line(aes(x = time, y = SMC_kill_rateout*1e4)) +
  #   facet_wrap(~arm+det) + 
  #   scale_y_log10() + 
  #   theme_classic()
  
  
  # Immunity over time by arm ####
  imm <- parasit %>%
    mutate(week = ceiling(time/7)) %>%
    group_by(arm, week) %>%
    summarise(innate_imm = 1 - mean(innate_imm),
              genadaptive_imm = 1 - mean(genadaptive_imm),
              varspecific_imm = 1 - mean(varspecific_imm),
              # time_imm = mean(st),
              growth = mean(growth),
              pb = mean(parasites))
  imm_agg <- parasit %>%
    group_by(arm) %>%
    summarise(innate_imm = mean(innate_imm),
              genadaptive_imm = mean(genadaptive_imm),
              varspecific_imm = mean(varspecific_imm),
              # time_imm = mean(st),
              growth = mean(growth)) %>%
    pivot_longer(innate_imm:growth, 
                 values_to = 'immunity_value',
                 names_to = 'immunity_type')
  # ggplot(imm_agg) + 
  #   geom_point(aes(x = immunity_type, y = immunity_value, group = arm, color = arm),
  #            position = 'dodge') +
  #   scale_color_manual(values = mycolors) +
  #   facet_wrap(~immunity_type, scales = 'free')
  
  inn <- ggplot(imm ) + 
    geom_point(aes(x = week, y = innate_imm, color = arm), alpha = 0.9) +
    geom_line(aes(x = week, y = innate_imm, color = arm), alpha=0.5) +
    scale_color_manual(values = mycolors) +
    theme_bw(base_size =  14) +#+ facet_wrap(~arm)
    labs(x = 'Weeks since follow-up',
         y = expression(paste('Innate immunity, 1 - S '[c],'(t)')))
  
  gen <- ggplot(imm ) + 
    geom_point(aes(x = week, y = genadaptive_imm, color = arm), alpha = 0.9) +
    geom_line(aes(x = week, y = genadaptive_imm, color = arm), alpha=0.5) +
    scale_color_manual(values = mycolors) +
    theme_bw(base_size =  14) +#+ facet_wrap(~arm)
    labs(x = 'Weeks since follow-up',
         y = expression(paste('General adaptive immunity, 1 - S '[m],'(t)')))
  
  var <- ggplot(imm ) + 
    geom_point(aes(x = week, y = varspecific_imm, color = arm), alpha = 0.9) +
    geom_line(aes(x = week, y = varspecific_imm, color = arm), alpha=0.5) +
    scale_color_manual(values = mycolors) +
    theme_bw(base_size =  14) +#+ facet_wrap(~arm)
    labs(x = 'Weeks since follow-up',
         y = expression(paste('var-specific immunity, 1 - S '[v],'(t)')))
  
  grow <- ggplot(imm ) + 
    geom_point(aes(x = week, y = growth, color = arm), alpha = 0.9) +
    geom_line(aes(x = week, y = growth, color = arm), alpha=0.5) +
    scale_color_manual(values = mycolors) +
    theme_bw(base_size =  14) + 
    labs(x = 'Weeks since follow-up',
         y = 'Growth rate (per two-day timestep)')
  ggsave(paste0('R:/Kelly/synergy_orderly/src/sim_cohort_generic/plots_infection_dynamics/growth_rate_by_arm_', type, outputsfolder, '.pdf'),
         plot = grow)
  
  timeimm <- ggplot(imm ) + 
    geom_point(aes(x = week, y = time_imm, color = arm)) +
    scale_color_manual(values = mycolors) +
    geom_line(aes(x = week, y = time_imm, color = arm), alpha=0.3) +
    theme_bw() #+ facet_wrap(~arm)
  
  
  library(cowplot)
  immunity_raw <- plot_grid(inn + theme(legend.position="none"), 
                            gen + theme(legend.position="none"),
                            # timeimm + theme(legend.position = 'none'),
                            var + theme(legend.position="none"), 
                            nrow = 1, 
                            labels = 'AUTO')
  # extract the legend from one of the plots
  legend <- get_legend(
    # create some space to the left of the legend
    inn + theme(legend.box.margin = margin(0, 0, 0, 12))
  )
  
  # add the legend to the row we made earlier. Give it one-third of 
  # the width of one plot (via rel_widths).
  immunity <- plot_grid(immunity_raw, legend, rel_widths = c(3, .3))
  # immunity
  ggsave(paste0('R:/Kelly/synergy_orderly/src/sim_cohort_generic/plots_infection_dynamics/immunity_by_arm_', type,  outputsfolder,'.pdf'),
         plot = immunity,
         height = 4, width = 13)
  
  
  
  
  # There is definitely less immunity in the combination intervention group which makes sense 
  
  # proportion of non-zero infections that become cases in RTSS vs no int group 
  # nonzero <- parasit %>%
  #   filter(mero_init_out != 0)
  # 
  # nonzero %>%
  #   janitor::tabyl(arm, det) %>%
  #   janitor::adorn_percentages()
  # basic: 
  # 88.5% of non-zero infections become cases in RTSS arm 
  # 98.3% of non-zero infections become cases in no int arm 
  
  # when innate immunity is 1, a higher prop (93.4%) of infections become cases in RTSS arm (also more in no int arm - 99.3%)
  # when genadaptive immunity is 1, av even higher prop (94.5%) of infections become cases in RTSS arm (lower in no int arm - 97.9%)
  
  # Plot incidence  ----
  inci_summarized <- inci %>%
    group_by(arm, yearmonth) %>%
    summarise(inci_median = median(incidence_per_1000pm, na.rm = TRUE),
              inci_lower = quantile(incidence_per_1000pm, 0.025, na.rm = TRUE),
              inci_upper = quantile(incidence_per_1000pm, 0.975, na.rm = TRUE))
  
  incidenceplot <- ggplot(inci_summarized %>% filter(yearmonth < '2018-06-01'))+
    geom_line(aes(x = as.Date(yearmonth), y = inci_median, color = arm), linewidth = 0.8) +
    geom_ribbon(aes(x = as.Date(yearmonth), ymin = inci_lower, ymax = inci_upper, fill = arm), alpha  = 0.4) +
    scale_x_date(breaks = '3 months',
                 labels = scales::label_date_short()) +
    scale_fill_manual(values = mycolors) +
    scale_color_manual(values =  mycolors,
                       breaks = c('both','none','rtss','smc','SMC delivery','RTS,S delivery'))+#c('#C44536','#772E25','#197278','#283D3B'))+
    labs(color = 'Intervention arm',
         x = 'Date',
         y = 'Incidence per 1000 person-months') +
    theme_bw(base_size = 14) + 
    facet_wrap(~factor(arm, levels = c('none','rtss','smc','both')))
  
  # Plot of the difference of model-predicted to expected by aggregation unit
  differenceplot <- inci_agg %>% 
    group_by(time_value) %>%
    summarise(difference_inci_averted_pred_exp = median(difference_inci_averted_pred_exp, na.rm = TRUE)) %>%
    ggplot(aes(x = as.Date(time_value), y = difference_inci_averted_pred_exp)) +
    # model estimated
    geom_point(size = 1) +
    geom_line() + 
    scale_x_date(date_breaks = '1 month',
                 labels = scales::label_date_short()) +
    # scale_color_manual(values = colors) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
    labs(
      x = NULL,#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Half-year',
      y = "Difference in model-predicted versus expected\ncases averted per 1000 people of\ncombination vs no intervention",
      shape = NULL, linetype = NULL,
      color = NULL#if(agg_unit == 'year') "Study year" else if (agg_unit == 'halfyear') 'Study half-year'
    ) +
    # scale_y_continuous(breaks = seq(-1,1,0.2)) + 
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = 'none')
  # when genadaptive immunity is always 1, there is less synergy (especially later in season)
  
  return(list(medianprbcs = medianprbcs,
              inn = inn, 
              gen = gen, 
              var = var, 
              grow = grow,
              immunity = immunity,
              incidenceplot = incidenceplot,
              differenceplot = differenceplot))
}

synergistic <- make_infection_dynamics_plots(outputsfolder = 'outputs_2026-03-12')
antagonistic <- make_infection_dynamics_plots(outputsfolder = 'outputs_2026-03-11_11' )
timingsyn <- make_infection_dynamics_plots(outputsfolder = 'outputs_2026-03-20')
timingant <- make_infection_dynamics_plots(outputsfolder = 'outputs_2026-03-20_2')
timingtest <- make_infection_dynamics_plots(outputsfolder = 'outputs_2026-03-10')

legend <- get_legend(
  # create some space to the left of the legend
  synergistic$medianprbcs + theme(legend.box.margin = margin(0, 0, 0, 12))
)

plot_grid(plot_grid(synergistic$medianprbcs + theme(legend.position = 'none'), 
          antagonistic$medianprbcs + theme(legend.position = 'none'), 
          nrow = 1,
          labels = 'AUTO'),
          legend,
          rel_widths = c(3, 0.3))
ggsave(paste0('R:/Kelly/synergy_orderly/src/sim_cohort_generic/plots_infection_dynamics/median_prbcs_combined', outputsfolder,'.pdf'),
       width = 10, height = 4.5)

plot_grid(plot_grid(synergistic$grow + theme(legend.position = 'none'), 
                    antagonistic$grow + theme(legend.position = 'none'), 
                    nrow = 1,
                    labels = 'AUTO'),
          legend,
          rel_widths = c(3, 0.3))
ggsave(paste0('R:/Kelly/synergy_orderly/src/sim_cohort_generic/plots_infection_dynamics/growth_rate_combined', outputsfolder,'.pdf'),
       width = 10, height = 4.5)


plot_grid(plot_grid(synergistic$inn + theme(legend.position = 'none'), 
                    synergistic$gen + theme(legend.position = 'none'), 
                    synergistic$var + theme(legend.position = 'none'), 
                    antagonistic$inn + theme(legend.position = 'none'), 
                    antagonistic$gen + theme(legend.position = 'none'),
                    antagonistic$var + theme(legend.position = 'none'),
                    nrow = 2,
                    labels = 'AUTO'),
          legend,
          rel_widths = c(3, 0.3))
ggsave(paste0('R:/Kelly/synergy_orderly/src/sim_cohort_generic/plots_infection_dynamics/immunity_combined', outputsfolder,'.pdf'),
       width = 13, height = 8.5)

