source("simulation_study_saved.R")

plot_eps <- 0.2

plot_df %>%
  filter(n >= 30, eps == plot_eps) %>%
  mutate(bdd = factor(bdd, c(T, F))) %>%
  ggplot(aes(x = n, y = RMSE, color = factor(mu), linetype = bdd)) +
  geom_line(linewidth = 0.75) + geom_point(size = 1) +
  labs(color = "Truth", linetype = "Constrained?",
       x = "Sample Size, n") +
  scale_x_log10() +
  scale_color_discrete(labels = c(expression("N(0.1, 0.04"^2*")"), 
                                  expression("N(0.5, 0.2"^2*")  "))) +
  theme_tufte() +
  theme(panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.5))

ggsave("Figures/Figure8.pdf", width = 6.75, height = 2.25, dpi=600, units = "in")



## COMPUTE ADDITIONAL QUANTITES REFERENCED IN SEC 4.1

# Linear Regression Slope of RMSE vs n
lm_RMSE <- plot_df %>%
  filter(eps == 0.2, !bdd) %>%
  lm(data = ., log(RMSE) ~ log(n) + mu)
lm_RMSE$coefficients[2] # Slope
summary(lm_RMSE)$r.squared # R-squared