# Analysis of summation range [f(t), t-4] over time
# where f(t) = lambda^(t-4) * (t-4) and lambda = 0.996

library(ggplot2)
library(dplyr)

# Parameters
lambda <- 0.996
t_max <- 100

# Create data frame
df <- data.frame(t = 1:t_max) %>%
  mutate(
    t_minus_4 = t - 4,
    f_t = ifelse(t <= 4, 1, lambda^(t-4) * (t-4)),
    range_start = ifelse(t <= 4, NA, f_t),
    range_end = ifelse(t <= 4, NA, t_minus_4),
    range_width = ifelse(t <= 4, 0, pmax(0, t_minus_4 - f_t)),
    sum_length = ifelse(t <= 4, 0, ceiling(range_width))  # Number of integer steps
  )

# Print key values at different time points
cat("Range [f(t), t-4] at key time points:\n")
cat("=====================================\n")
key_times <- c(1,2,3,4,5,6,7,8,9, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100)
for(t_val in key_times) {
  row <- df[df$t == t_val, ]
  cat(sprintf("t=%d: f(t)=%.3f, t-4=%d, range=[%.3f, %d], width=%.3f\n", 
              t_val, row$f_t, row$t_minus_4, row$range_start, row$range_end, row$range_width))
}

# Plot 1: f(t) and t-4 over time
p1 <- ggplot(df %>% filter(t > 4), aes(x = t)) +
  geom_line(aes(y = f_t, color = "f(t) = λ^(t-4)*(t-4)"), size = 1) +
  geom_line(aes(y = t_minus_4, color = "t-4"), size = 1) +
  labs(title = "Summation Bounds Over Time",
       subtitle = "λ = 0.996",
       x = "Time t", y = "Value",
       color = "Function") +
  theme_minimal() +
  scale_color_manual(values = c("f(t) = λ^(t-4)*(t-4)" = "blue", "t-4" = "red"))

# Plot 2: Range width over time
p2 <- ggplot(df %>% filter(t > 4), aes(x = t, y = range_width)) +
  geom_line(color = "#3EC74F", size = 1) +
  geom_point(data = df %>% filter(t %in% key_times), 
             aes(x = t, y = range_width), color = "#1D5C25", size = 2) +
  scale_x_continuous(breaks = seq(0,100, 10), expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = seq(0,40, 2), expand = c(0.01,0.01))+
  labs(#title = "var-specific Range Width Over Time",
       # subtitle = "Width = (t-4) - f(t)",
       x = "Timestep", y = "Window width (2-day timesteps)") +
  theme_minimal()
ggsave(filename = 'R:/Kelly/synergy_orderly/figures/f_t_versus_t-4.pdf', p2, height = 4, width = 6)

# Plot 3: Visualization of the range as filled area
p3 <- ggplot(df %>% filter(t > 4 & t <= 50), aes(x = t)) +
  geom_ribbon(aes(ymin = f_t, ymax = t_minus_4), alpha = 0.3, fill = "green") +
  geom_line(aes(y = f_t, color = "f(t)"), size = 1) +
  geom_line(aes(y = t_minus_4, color = "t-4"), size = 1) +
  labs(title = "Summation Range Visualization (t = 5 to 50)",
       subtitle = "Shaded area shows range [f(t), t-4]",
       x = "Time t", y = "Value",
       color = "Boundary") +
  theme_minimal() +
  scale_color_manual(values = c("f(t)" = "blue", "t-4" = "red"))

p1
p2
p3



# Summary statistics
summary_stats <- df %>% 
  filter(t > 4) %>%
  summarise(
    min_width = min(range_width, na.rm = TRUE),
    max_width = max(range_width, na.rm = TRUE),
    final_width = last(range_width),
    final_f_t = last(f_t),
    final_t_minus_4 = last(t_minus_4)
  )
summary_stats

# Show how many discrete time steps are included in the sum at different times
# Number of discrete steps in summation:
for(t_val in key_times) {
  row <- df[df$t == t_val, ]
  start_step <- ceiling(row$range_start)
  end_step <- row$range_end
  n_steps <- max(0, end_step - start_step + 1)
  cat(sprintf("t=%d: steps from %d to %d = %d steps total\n", 
              t_val, start_step, end_step, n_steps))
}
