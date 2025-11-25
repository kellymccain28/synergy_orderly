# Script for visualizing/prepping fitting of RTSS parameters 
# simulation script is fit_rtss.R

source(paste0(path, '/shared/rtss.R'))

# Prepare the "observed" data to which we will compare the cohort sim model output
# Efficacy against infection from White 2015
nsim <- 100
ts <- seq(1,365*3)
phases <- ifelse(ts < 366, 1,
                 ifelse(ts < 729, 2, 3))
csp <- list()
for(i in 1:nsim){
  csp[[i]] <- antibody_titre(t = ts,
                     phase = phases,
                     peak1 = c(621,0.35) ,
                     peak2 = c(277, 0.35),
                     peak3 = c(277,0.35),
                     duration1 = c(45,16),
                     duration2 = c(591,245),
                     rho1 = c(2.37832, 1.00813),
                     rho2 = c(1.034, 1.027),
                     rho3 = c(1.034, 1.027))
}
csp_df <- data.frame(do.call(rbind, csp))
colnames(csp_df) <- ts
csp_df$sim <- 1:nsim

csp_long <- pivot_longer(csp_df, cols = -sim,
                         names_to = "time", values_to = "titre")
csp_long$time <- as.numeric(csp_long$time)

# Plot
ggplot(csp_long, aes(x = time, y = titre, group = sim)) +
  geom_line(alpha = 0.3) +
  theme_minimal() + scale_y_log10() +
  labs(x = "Time", y = "Antibody Titre")

ve <- lapply(csp, vaccine_eff)
ve_df <- data.frame(do.call(rbind, ve))
colnames(ve_df) <- ts
ve_df$sim <- 1:nsim

ve_long <- pivot_longer(ve_df, cols = -sim,
                         names_to = "time", values_to = "ve")
ve_long$time <- as.numeric(ve_long$time)
ve_long <- ve_long %>% ungroup() %>%
  mutate(weeks_since_rtss = ceiling(time/7)) %>%
  # filter(weeks_since_rtss < 52) %>%
  group_by(weeks_since_rtss, sim) %>%
  summarize(ve_inf = mean(ve))

# Get summary of all efficacy values over time 
ve_long <- ve_long %>% 
  group_by(weeks_since_rtss) %>%
  mutate(observed_efficacy = median(ve_inf),
            efficacy_lower = quantile(ve_inf, 0.025),
            efficacy_higher = quantile(ve_inf, 0.975))

# Plot
ggplot(ve_long) +
  geom_line(aes(x = weeks_since_rtss, y = ve_inf, group = sim), color = 'lightblue') +
  geom_line(aes(x = weeks_since_rtss, y = observed_efficacy), color = 'darkgreen', linewidth = 2) + 
  # geom_errorbar(aes(x = weeks_since_rtss, ymin = efficacy_lower, ymax = efficacy_higher)) +
  theme_minimal() +
  labs(x = "Weeks since RTSS", y = "Efficacy")

saveRDS(ve_long %>% select(weeks_since_rtss, observed_efficacy) %>% distinct(), file = paste0(path, '/src/fit_rtss/observed_rtss_efficacy.rds'))
