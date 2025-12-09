library(ggplot2)
library(umbrella)

# Read in the generic seasonal curve of the probability of an infectious bite 
prob_bite_generic <- readRDS("R:/Kelly/synergy_orderly/archive/fit_rainfall/20251203-121008-03d1b614/prob_bite_generic.rds")

ggplot(prob_bite_generic %>% filter(year == 2018)) +
  geom_line(aes(x = date, y = prob_infectious_bite, color = '#AE0939'), linewidth = 2)

seasonal <- c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919)#https://github.com/mrc-ide/mlgts/blob/master/data/seasonal.rda
seasonal <- set_names(seasonal, c('g0','g1','g2','g3','h1','h2','h3'))
generic_seas <- fourier_predict(coef = seasonal, t = 1:365, floor = 0.001)
generic_seas <- generic_seas %>%
  mutate(date = as.Date(t, origin = '2017-01-01'),
         date_clin = date + 45) # to approximate the lag between rainfall and probability of a bite 

smc_lines <- data.frame(
  xint = c(as.Date('2017-04-01') + 115 + c(0,30,60,90)),
  color = '#4D9DE0'
)
smc_shaded <- data.frame(
  xint = c(as.Date('2017-04-01') + 115 + c(0,30,60,90, 120)),
  color = '#4D9DE0'
)
# metadata_df$vaccination_day[1] = 90
rtss_lines <- data.frame(
  xint = c(as.Date('2017-04-01') + c(75,75-30, 75-60)),
  color = '#59114D'
)

ggplot(generic_seas) + 
  # geom_line(data = prob_bite_generic %>% filter(year == 2017), aes(x = date, y = prob_infectious_bite*5, color = '#AE0939'), linewidth = 2)+
  geom_line(aes(x = date_clin, y = profile, color = "Probability of\ninfectious bite"), linewidth = 2) +
  geom_vline(data = rtss_lines, aes(xintercept = xint, color = 'Vaccine doses'), linetype = 2) +
  geom_vline(data = smc_lines, aes(xintercept = xint, color = 'SMC rounds'), linetype = 2) +
  geom_area(data = smc_shaded, aes(x = xint, y = 1), fill = colorspace::lighten('#4D9DE0', amount = 0.5), alpha = 0.15) +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_blank(),
        legend.position = c(0.9,0.8)) + 
  scale_x_date(breaks = '1 month',
               labels = scales::label_date_short())+
  scale_color_manual(values =  c('SMC rounds' = '#4D9DE0',
                                 'Vaccine doses' = '#470024',
                                 "Probability of\ninfectious bite" = '#AE0939'))+
  labs(y = NULL,
       x = NULL,
       color = NULL) 

ggsave('R:/Kelly/synergy_orderly/figures/generic_seasonal_interventions.pdf')
 
# 
ggplot(outputs_generic %>%
         filter(date < as.Date('2018-01-01'))) +
  geom_point(aes(x = date, y = n_bitten)) +
  geom_line(data = generic_seas, aes(x = date, y = profile*1000)) +
  scale_x_date(breaks = '1 month') + 
  theme(axis.text.x = element_text(angle = 90))
