# script to quickly look at dynamics of smc infection over time 
observed_efficacy <- read.csv("R:/Kelly/synergy_orderly/shared/smc_fits_hayley.csv")

# cohort runs using 'fitted' parameters for smc 
test_fitted_params_smc <- readRDS("R:/Kelly/synergy_orderly/src/fit_smc/outputs/test_fitted_params_smc_0611.rds")
paras <- test_fitted_params_smc[[1]]$parasitemia

p1 <- paras %>%
  filter(rid == 7)
ggplot(p1) + 
  geom_point(aes(x = time_ext, y = SMC_kill_rateout, color = as.factor(day1_BSinfection), shape = as.factor(det))) 
ggplot(p1) + 
  geom_point(aes(x = time_ext, y = prob_smckill, color = as.factor(day1_BSinfection), shape = as.factor(det))) 
ggplot(p1) + 
  geom_point(aes(x = time_ext, y = numkillSMC, color = as.factor(day1_BSinfection)), shape = 4) +
  geom_point(aes(x = time_ext, y = parasites, color = as.factor(day1_BSinfection), shape = as.factor(det)))+
  scale_y_log10()

table(p1$det)/100

ggplot(p1%>% filter(day1_BSinfection == 65)) + 
  geom_point(aes(x = time, y = prob_smckill, color = as.factor(day1_BSinfection)), shape = 4) 

ggplot(p1%>% filter(day1_BSinfection == 65)) + 
  geom_point(aes(x = time, y = parasites, color = as.factor(day1_BSinfection), shape = as.factor(det)))+
  scale_y_log10()


# Test SMC with different times of infection relative to SMC 
intervalsmc_days = 30 
smc_dose_days <- 0#c(seq(runpars$season_start_day, runpars$season_start_day + 120 - 1, intervalsmc_days),
# seq(runpars$season_start_day + 365, runpars$season_start_day + 365 + 120 - 1, intervalsmc_days),
# seq(runpars$season_start_day + 365*2, runpars$season_start_day  + 365*2 + 120 - 1, intervalsmc_days))
max_SMC_kill_rate=3.07#2.43#5#6
lambda = 13#11.6#28#11.6
kappa = 0.454#0.173#5.42#1.81
ts = 100
ratetest <- max_SMC_kill_rate * exp(-(seq(1:ts)/lambda)^kappa)
plot(ratetest)
plot(1-exp(-ratetest * 1))

threshold = 5000
n_particles = 100

get_smc_outputs <- function(inf_start){
  tt <- seq(0:(365-inf_start))
  
  time_since_smc <- lapply(list(smc_dose_days), 
                           calc_time_since_dose, 
                           days = 0:(ts*2 -1))
  smckillvec <- lapply(time_since_smc,
                       function(.){
                         kill <- max_SMC_kill_rate * exp(-(./ lambda)^kappa)  # calculate kill rate with hill function
                         kill[is.na(kill)] <- 0                                       # change NAs (when there is no SMC) to 0 
                         # kill <- kill[seq_along(kill) %% 2 == 1]
                         # Sum every two consecutive days 
                         n <- length(kill)
                         if (n %% 2 == 1) kill <- c(kill, kill[n])                         # pad with last value if odd length
                         colSums(matrix(kill, nrow = 2))
                       })
  smc_killvec <- unlist(smckillvec)
  smckillvec_subset <- smc_killvec[floor((inf_start) / 2) :length(smc_killvec)]
  smckillvec_subset <- list(c(smckillvec_subset, rep(0, ts-length(smckillvec_subset))))
  
  df <- run_model(n_particles = n_particles,
                   n_threads = 4L,
                   PEV_on = 0,
                   SMC_on = 1,
                   tt= tt,
                   det_mode = FALSE,
                   t_inf_vax = 0,
                   VB = 1e6,
                   SMC_time = seq(0, length(smckillvec[[1]])-1,1),
                   SMC_kill_vec = smckillvec_subset,
                  infection_start_day = inf_start
  ) %>%
    format_data(tt= tt,
                infection_start_day = inf_start,
                n_particles = n_particles) %>%
    make_plots()
  
  dd <- df[[6]] %>% ungroup() %>%select(run, detectable) %>%
    distinct() %>% janitor::tabyl(detectable) %>%
    mutate(infection_start_day = inf_start)
  return(list('efficacy' = dd, 
              'df'= df[[6]] %>% mutate(arm = 'smc', inf_start = inf_start)))
}
get_none_outputs <- function(inf_start){
  tt <- seq(0:(365-inf_start))
  
  df <- run_model(n_particles = n_particles,
                  n_threads = 4L,
                  PEV_on = 0,
                  SMC_on = 0,
                  tt= tt,
                  det_mode = FALSE,
                  t_inf_vax = 0,
                  VB = 1e6,
                  SMC_time = seq(0,max(tt),1),
                  SMC_kill_vec = rep(0,max(tt)+1),
                  infection_start_day = inf_start
  ) %>%
    format_data(tt= tt,
                infection_start_day = inf_start,
                n_particles = n_particles) %>%
    make_plots()
  
  dd <- df[[6]] %>% ungroup() %>% select(run, detectable) %>%
    distinct() %>% janitor::tabyl(detectable) %>%
    mutate(infection_start_day = inf_start)
  
  return(list('efficacy' = dd, 
              'df'= df[[6]] %>% mutate(arm = 'smc', inf_start = inf_start)))
}

test_none <- lapply(1:70, get_none_outputs)
eff_none <- purrr::map_df(test_none, 'efficacy')
df_none <- purrr::map_df(test_none, 'df')
nonon <- bind_rows(eff_none) %>%
  filter(detectable == 'not detectable')%>%
  rename(n_none = n,
         perc_none = percent) %>%
  mutate(n_inf_none = 100-n_none)

test_div2 <- lapply(1:70, get_smc_outputs)
eff_smc <- purrr::map_df(test_div2, 'efficacy')
df_smc <- purrr::map_df(test_div2, 'df')
fdaf <- eff_smc %>% 
  filter(detectable == 'not detectable') %>%
  mutate(n_inf = 100 - n)

lalal <- nonon %>% left_join(fdaf, by = c('infection_start_day', 'detectable'))
difff <- lalal %>% 
  mutate(eff = 1-(n_inf / n_inf_none))
ggplot(difff)+
  geom_point(aes(x = 1:70, y = eff))+
  geom_line(data = observed_efficacy, aes(x = day_since_smc, y = efficacy), 
            color = 'blue', linewidth = 1, linetype = 1) +
  # geom_line(aes(x = 1:70, y = ratetest/max(ratetest))) +
  ylim(c(0,1)) + 
  theme_bw() + 
  labs(caption = str_glue(max_SMC_kill_rate, ', ', lambda,', ',kappa))

ggplot(difff) + 
  geom_line(aes(x = 1:70, y = n_inf_none), color = 'black') + 
  geom_line(aes(x = 1:70, y = n_inf), color = 'darkred')

df_smc %>%
  ggplot() + 
  geom_line(aes(x = time, y = smcrate, group = inf_start)) + 
  xlim(c(0,70))
