library(ggplot2)
library(umbrella)


# Read in the generic seasonal curve 
prob_bite_generic <- readRDS("R:/Kelly/synergy_orderly/archive/fit_rainfall/20251120-190203-3e232d5d/prob_bite_generic.rds")

# ggplot(prob_bite_generic %>% filter(year == 2020)) + 
#   geom_line(aes(x = date, y = prob_infectious_bite, color = '#AE0939'), linewidth = 2)

seasonal <- c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919)
seasonal <- set_names(seasonal, c('g0','g1','g2','g3','h1','h2','h3'))
generic_seas <- fourier_predict(coef = seasonal, t = 1:365, floor = 0.001)
generic_seas <- generic_seas %>%
  mutate(date = as.Date(t, origin = '2020-01-01'))
ggplot(generic_seas) + 
  geom_line(aes(x = date, y = profile), color = '#AE0939', linewidth = 2) +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_blank()) + 
  scale_x_date(breaks = '1 month',
               labels = scales::label_date_short())+
  labs(y = NULL,
       x = NULL)
ggsave('R:/Kelly/synergy_orderly/figures/generic_seasonal.pdf')
