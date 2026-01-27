
library(dplyr)

recrudescence <- p %>% 
  filter(rid == 112 & day1_BSinfection == 6) %>%
  group_by(rid, day1_BSinfection) %>%
  arrange(rid, day1_BSinfection, time_withinhost) %>%
  mutate(
    # Identify when parasites are above/below 5000
    above_5000 = parasites > 5000,
    below_5000 = parasites <= 5000,
    # Detect all potential crossings from below to above
    crossed_up = above_5000 & lag(below_5000, default = FALSE)
  ) %>%
  summarise(
    threshold_day = first(threshold_day),
    # Get all potential crossing days
    all_crossings = list(time_withinhost2[crossed_up]),
    # Filter crossings: first must be 7+ days after threshold, 
    # subsequent must be 7+ days after previous
    valid_crossings = list({
      crossings <- all_crossings[[1]]
      if (length(crossings) == 0) return(numeric(0))
      
      valid <- numeric()
      last_valid <- threshold_day[1]
      
      for (day in crossings) {
        if (day >= last_valid + 7) {
          valid <- c(valid, day)
          last_valid <- day
        }
      }
      valid
    }),
    second_infection_day = ifelse(length(valid_crossings[[1]]) >= 1, 
                                  valid_crossings[[1]][1], 
                                  NA),
    third_infection_day = ifelse(length(valid_crossings[[1]]) >= 2, 
                                 valid_crossings[[1]][2], 
                                 NA),
    .groups = "drop"
  )
# View results
recrudescence %>%
  filter(!is.na(second_infection_day))

# Count how many have second/third infections
cat("Second infections:", sum(!is.na(recrudescence$second_infection_day)), "\n")
cat("Third infections:", sum(!is.na(recrudescence$third_infection_day)), "\n")

p5 <- p %>% filter(rid == 5 & day1_BSinfection == 767)
ggplot(p5) + geom_line(aes (x = time_withinhost, y = parasites)) + 
  geom_hline(yintercept = 5000, linetype = 2) +
  scale_y_log10()+xlim(c(0,75))

p112 <- p %>% filter(rid == 112 & day1_BSinfection == 6)
ggplot(p112) + 
  geom_line(aes (x = time_withinhost2, y = parasites)) + 
  geom_hline(yintercept = 5000, linetype = 2) + 
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0,130,7), limits = c(0,120))+
  geom_vline(xintercept = unlist(recrudescence$valid_crossings), linetype = 4) + 
  theme_bw()
