# script to make figures for thesis 

source('shared/rtss.R')
# Distribution of sporozoites and merozoites with and without vaccination ----

# Get antibodies 
ts <- seq(0,365*3)
tboost1 = 365
tboost2 = 730
phases <- ifelse(ts < tboost1, 1,
                 ifelse(ts < tboost2, 2, 3))
# values from hogan : https://github.com/mrc-ide/rtss_vacc_antibody_model/blob/main/main.R
csp <- antibody_titre(t = ts,
                     phase = phases,
                     peak1 = c(621,0.35) ,
                     peak2 = c(277, 0.35),
                     peak3 = c(277,0.35),
                     duration1 = c(45,16),
                     duration2 = c(591,245),
                     rho1 = c(2.37832, 1.00813),
                     rho2 = c(1.034, 1.027),
                     rho3 = c(1.034, 1.027),
                     t_boost1 = tboost1,
                     t_boost2 = tboost2)
#Dose response
ab <- 0#
beta_ab <- 6.62 #(1.34-16.29) # anti-CSP titre for 50% reduction in spz survival prob microgram/mL
alpha_ab <- 1.32 #(0.85-1.77) # shape parameter for antibody dose-response

DR <- 1 / (1 + (ab / beta_ab)^alpha_ab) # prob of survival of single spz
# ggplot()+geom_point(aes(x = ts, y = ab))+scale_y_log10()
plot(DR)

# Parameters for spz model initial merozoites 
# estimated and fixed parameters  
n <- 150 #n, mean number of successful spz per challenge / 5 bites; neg bin
sigma_n <- 194 #, sigman sd of number of successful spz per challenge / 5 bites
mu <- 2136 #30000 #10.1371/journal.pcbi.1005255 as assumed by Hayley #2136 # mean number of merozoites released per sporozoite in Michael's model; gamma distributed
sigma_mu <- 4460 #71427 #4460 from Michael's model # sd of number of merozoites released per sporozoite; gamma distributed

# Parameters for Negative Binomial distribution
r <- n^2 / (sigma_n^2 - n)
p <- r / (n*DR + r) #n*DR / (n*DR + r)- this is 1- r / (n*DR + r) # should be probability of success aka prob of survival 

# Draw number of successful sporozoites 
xspz <- seq(0, 100)
kspz <- dnbinom(xspz, r / 5, p)
p1 <- ggplot()+
  geom_area(aes(x = xspz, y = kspz), 
            color = '#35409C', alpha = 0.4, fill = '#636DCA', linewidth = 1) + 
  geom_vline(aes(xintercept = n), color = '#35409C', linetype = 2) +
  scale_x_continuous(expand = c(0,5)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_minimal(base_size = 14) + 
  theme(axis.text.y = element_blank()) + 
  labs(x = 'Number of sporozoites',
       y = 'Density')
p1

# Merozoites per sporozoite (k = 1)
# Parameters for Gamma distribution
kspz2=1
theta <- sigma_mu^2 / mu  # Scale parameter (theta)
gamma_shape <- (kspz2 * mu) / theta
# find total merozoites for from k sporozoites
mero_init <- rgamma(10000, shape = gamma_shape, scale = theta) 
# ggplot()+geom_density(aes(x = mero_init))
mean(mero_init)

# or 
x_mero <- seq(0, 10000, length.out = 1000)
mero_init2 <- dgamma(x_mero, shape = gamma_shape, scale = theta)

p2 <- ggplot()+
  geom_area(aes(x = x_mero, y = mero_init2), 
            color = '#9C2007', alpha = 0.4, fill = '#FAA18F', linewidth = 1)+ 
  geom_vline(aes(xintercept = mu), color = '#9C2007', linetype = 2) +
  scale_x_continuous(expand = c(0,200)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_minimal(base_size = 14) + 
  theme(axis.text.y = element_blank()) + 
  labs(x = 'Number of merozoites per sporozoite',
       y = 'Density')
p2

# Total merozoites for a range of sporozoite values 
# scale the gamma mean and sd for the number of sporozoites (mean n)
n = 15
mean_total <- n * mean(mero_init) # Mean number of spz * mean number of mero per spz
sd_total <- sqrt(n) * sd(mero_init)

scale_total <- sd_total^2 / mean_total
shape_total <- mean_total / scale_total

mero150spz <- dgamma(seq(1e4,2e5, length.out= 1e3), shape = shape_total, scale = scale_total)

n2 = 30
mean_total2 <- n2 * mean(mero_init) # Mean number of spz * mean number of mero per spz
sd_total2 <- sqrt(n2) * sd(mero_init)

scale_total2 <- sd_total2^2 / mean_total2
shape_total2 <- mean_total2 / scale_total2

mero30spz <- dgamma(seq(1e4, 2e5, length.out= 1e3), shape = shape_total2, scale = scale_total2)

p3 <- ggplot() + 
  geom_area(aes(x = seq(1e4, 2e5, length.out = 1e3), y = mero30spz, color = '15 sporozoites'),
            # color = '#7D359C',
            alpha = 0.4, fill = '#BD85D5', linewidth = 1)+
  geom_area(aes(x = seq(1e4, 2e5, length.out = 1e3), y = mero150spz, color = '30 sporozoites'),
            # color = '#1F843F', 
            alpha = 0.4, fill = '#28A951', linewidth = 1)+
  geom_vline(aes(xintercept = mean_total2), color = '#7D359C', linetype = 2) +
  geom_vline(aes(xintercept = mean_total), color = '#1F843F', linetype = 2) +
  scale_x_continuous(expand = c(0,6000),
                     labels = scales::scientific_format(digits=1)) +
  scale_y_continuous(expand = c(0,0.0000003)) +
  scale_color_manual(values = c('15 sporozoites' = '#BD85D5',
                                '30 sporozoites' = '#28A951'),
                     name = 'k sporozoites\ninitiating infection') +
  theme_minimal(base_size = 14) + 
  theme(axis.text.y = element_blank(),
        legend.position = c(0.82, 0.85)) + 
  labs(x = 'Total number of merozoites',
       y = 'Density')
p3

all_plots <- cowplot::plot_grid(p1, p2, p3, nrow = 1, 
                                labels = c('A','B','C'))
ggsave(filename = 'R:/Kelly/synergy_orderly/figures/spz_mero_distributions.pdf', all_plots,
       height = 4.5, width = 13)


# Antibodies versus vaccine effectiveness


# ve <- vaccine_eff(ab)



# 