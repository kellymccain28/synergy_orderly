


all_results <- rbind(tidy_results %>% mutate(type = 'trial'),
                     tidy_results_model %>% mutate(type = 'model')) %>%
  filter(term %in% c("Both vs. RTSS",'RTSS vs. SMC','Both vs. SMC'))

ggplot(all_results)+
  geom_point(aes(x = term, y = VE, group = type, color = type), #color = 'darkgreen', 
             position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(x = term, ymin = VE_lower, ymax = VE_upper, group = type, color = type), 
                width = 0.2, #color = 'darkgreen', 
                position = position_dodge(width = 0.3),
                linewidth = 1) +
  scale_color_manual(values = c('darkgreen','orange')) +
  geom_hline(aes(yintercept = 0)) +
  labs(y = "Efficacy",
       x = NULL,
       color = NULL,
       caption = 'Green is model, red is trial') +
  scale_y_continuous(breaks = seq(-0.2, 1, 0.1),
                     limits = c(-0.2,1),
                     labels = scales::percent) + 
  theme_bw(base_size = 16)
