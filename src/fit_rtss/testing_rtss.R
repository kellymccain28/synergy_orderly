# Script to test out different values of alpha and beta for the rtss calibration 
# 7.1.2026

source('shared/rtss.R')

ts <- seq(0,365*2)
# phases <- ifelse(ts < 365, 1,
#                  ifelse(ts < 730, 2, 3))
# ab <- antibody_titre(t = ts,
#                      phase = phases,
#                      peak1 = c(621,0.35) ,
#                      peak2 = c(277, 0.35),
#                      peak3 = c(277,0.35),
#                      duration1 = c(45,16),
#                      duration2 = c(591,245),
#                      rho1 = c(2.37832, 1.00813),
#                      rho2 = c(1.034, 1.027),
#                      rho3 = c(1.034, 1.027),
#                      t_boost1 = 365,
#                      t_boost2 = 730)
ab <- antibody_titre_det(t = ts, 
                   t_boost1 = 365, 
                   t_boost2 = 730)
plot(ab)

# dose response 
spzeff1 <- p_spz_surv(ab, beta_ab = 5.83, alpha_ab = 2)
spzeff2 <- p_spz_surv(ab, beta_ab = 5.83, alpha_ab = 1.32)
spzeff3 <- p_spz_surv(ab, beta_ab = 5.83, alpha_ab = 1.6)

spzeff4 <- p_spz_surv(ab, beta_ab = 5, alpha_ab = 1.6)
spzeff5 <- p_spz_surv(ab, beta_ab = 5, alpha_ab = 2)
spzeff6 <- p_spz_surv(ab, beta_ab = 5, alpha_ab = 1.32)

DR <- p_spz_surv(ab, beta_ab = 6.62, alpha_ab = 1.32)
spzeff8 <- p_spz_surv(ab, beta_ab = 6.62, alpha_ab = 2)
spzeff9 <- p_spz_surv(ab, beta_ab = 6.62, alpha_ab = 1.6)

spzeff10 <- p_spz_surv(ab, beta_ab = 5.83, alpha_ab = 1.48)
spzeff10 <- p_spz_surv(ab, beta_ab = 2.7, alpha_ab = 1.8)

DR <- p_spz_surv(ab, beta_ab = 2.46, alpha_ab = 1.59)

n <- 150 #n, mean number of successful spz per challenge ; neg bin
sigma_n <- 194 
r <- (n^2) / (sigma_n^2 - n) 
p <-  n*DR / (n*DR + r)#vmin + (1 - vmin) *  
n_dr = n*DR
plot(p)
plot(n_dr)

plot(dnbinom(1:400, size = r, mu = n_dr))

ggplot() + 
  # geom_line(aes(x = ts, y = 1 - spzeff1, color = '5.83, 2'), linewidth = 1) +
  geom_line(aes(x = ts, y = 1 - spzeff2, color = '5.83, 1.38'), linewidth = 1) +
  # geom_line(aes(x = ts, y = 1 - spzeff3, color = '5.83, 1.6'), linewidth = 1) +
  # geom_line(aes(x = ts, y = 1 - spzeff4, color = '5, 1.6'), linewidth = 1) +
  # geom_line(aes(x = ts, y = 1 - spzeff5, color = '5, 2'), linewidth = 1) +
  # geom_line(aes(x = ts, y = 1 - spzeff6, color = '5, 1.32'), linewidth = 1) +
  geom_line(aes(x = ts, y = 1 - spzeff7, color = '6.62, 1.32'), linewidth = 1) +
  # geom_line(aes(x = ts, y = 1 - spzeff8, color = '6.62, 2'), linewidth = 1) +
  # geom_line(aes(x = ts, y = 1 - spzeff9, color = '6.62, 1.6'), linewidth = 1) +
  geom_line(aes(x = ts, y = 1 - spzeff10, color = '10'), linewidth = 1) +
  geom_line(aes(x = ts, y = 1 - spzeff11, color = 'new'), linewidth = 1) +
  theme_bw()

# higher beta steepens
# higher alpha flattens it and increases starting value 
# combo of higher alpha and lower beta flattens the curve 

beta_ab <- seq(4, 6, by = 0.2)
alpha_ab <- seq(1.2, 2, by = 0.1)
parvals <- crossing(beta_ab, alpha_ab)

pars <- params_df %>% #parvals %>%
  rowwise() %>%
  mutate(dr = list(p_spz_surv(ab, beta_ab, alpha_ab))) %>%
  ungroup() %>%
  mutate(ts = list(ts)) %>%  # Add ab_vector to each row
  unnest(cols = c(ts, dr)) 

# pars <- pars %>% mutate(group_id = rep(1:nrow(pars), each = length(ab))) 

ggplot(pars %>% filter(alpha_ab > 1.2), aes(x = ts, y = 1 - dr, 
                                            color = alpha_ab)) +
  geom_point() +
  geom_line(data = pars %>% filter(beta_ab == 6.62), aes(x = ts, y = 1-dr), color = 'red') +
  geom_line(data = pars %>% filter(beta_ab == 5.83), aes(x = ts, y = 1-dr), color = 'red') +
  scale_color_viridis_c() +  # or scale_color_gradient() for custom colors
  theme_minimal() + 
  ylim(c(0.9, 1))






# Parameters from the model
n <- 150           # mean number of surviving sporozoites (unvaccinated)
sig_n <- 194       # standard deviation
r <- (n^2)/(sig_n^2 - n) /5 # shape parameter

# above is the same as :
n <- 150/5
sig_n <- 194/sqrt(5)
r <- (n^2)/(sig_n^2 - n)


# Range of possible values for k (number of surviving sporozoites)
k <- 0:500

# Calculate PMF for unvaccinated (baseline)
DR = 1
mu_unvax <- n 
p <- n*DR / (n*DR + r) 
pmf_unvax <- dnbinom(k, size = r, mu = mu_unvax)
pmf_unvax <- dnbinom(k, size = r, prob = 1-p)

# Calculate PMF for vaccinated scenarios with different antibody levels
# Example: DR = 0.5 (moderate vaccine effect)
DR_moderate <- 0.5
mu_moderate <- n * DR_moderate
p <- n*DR / (n*DR + r) 
pmf_vax_moderate <- dnbinom(k, size = r, mu = mu_moderate)
pmf_vax_moderate <- dnbinom(k, size = r, prob = 1-p)

# Example: DR = 0.2 (strong vaccine effect)
DR_strong <- 0.2
mu_strong <- n * DR_strong
pmf_vax_strong <- dnbinom(k, size = r, mu = mu_strong)

# Plot
par(mfrow = c(1, 1))
plot(k, pmf_unvax, type = "h", lwd = 2, col = "black",
     xlab = "Number of surviving sporozoites (k)",
     ylab = "Probability P(X = k)",
     main = "PMF of Sporozoite Survival (mu parameterization)",
     ylim = c(0, max(pmf_unvax)))
lines(k, pmf_vax_moderate, type = "h", lwd = 2, col = "blue")
lines(k, pmf_vax_strong, type = "h", lwd = 2, col = "red")
legend("topright", 
       legend = c(paste0("Unvaccinated (μ=", mu_unvax, ")"),
                  paste0("Moderate vaccine (μ=", mu_moderate, ")"),
                  paste0("Strong vaccine (μ=", mu_strong, ")")),
       col = c("black", "blue", "red"), 
       lwd = 2)

# Add vertical lines for means
abline(v = mu_unvax, col = "black", lty = 2)
abline(v = mu_moderate, col = "blue", lty = 2)
abline(v = mu_strong, col = "red", lty = 2)
