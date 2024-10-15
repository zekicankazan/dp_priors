source("simulation_study_saved.R")

plot_eps <- 0.2

# Create Top Panel
p1 <- plot_df %>%
  filter(n >= 30, eps == plot_eps) %>%
  mutate(bdd = factor(bdd, c(T, F))) %>%
  ggplot(aes(x = n, y = CI_len, color = factor(mu), linetype = bdd)) +
  geom_line(linewidth = 0.75) + geom_point(size = 1) +
  labs(color = "Truth", linetype = "Constrained?", x = "Sample Size, n", y = "Avg. Interval Length") +
  scale_x_log10() +
  scale_color_discrete(labels = c(expression("N(0.1, 0.04"^2*")"), expression("N(0.5, 0.2"^2*")  "))) +
  theme_tufte()  +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

# Create Bottom Panel
p2 <- plot_df %>%
  filter(n >= 30, eps == plot_eps) %>%
  mutate(bdd = factor(bdd, c(T, F))) %>%
  ggplot(aes(x = n, y = coverage, color = factor(mu), linetype = bdd)) +
  geom_line(linewidth = 0.75) + geom_point(size = 1) +
  labs(color = "Truth", linetype = "Constrained?", x = "Sample Size, n", y = "Coverage") +
  scale_x_log10() +
  scale_color_discrete(labels = c(expression("N(0.1, 0.04"^2*")"), expression("N(0.5, 0.2"^2*")  "))) +
  theme_tufte()

# Create Figure 2
p1/p2 + plot_layout(guides = "collect") & 
  theme(plot.margin=grid::unit(c(1,0.1,4,0), "mm"), legend.position="bottom",
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(-5,0,-25,0),
        legend.box="vertical", legend.spacing.y = unit(0, "mm"),
        panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.5))

ggsave("Figures/Figure2.pdf", width = 3, height = 3.6, dpi=600, units = "in")


## COMPUTE ADDITIONAL QUANTITES REFERENCED IN SEC 4.1

# Linear Regression Slope of CI_length vs n
lm_CI <- plot_df %>%
  filter(eps == 0.2, !bdd) %>%
  lm(data = ., log(CI_len) ~ log(n) + mu)
lm_CI$coefficients[2] # Slope
summary(lm_CI)$r.squared # R-squared

# Max Interval Length For Constrained Analysis
plot_df %>%
  filter(bdd, eps == 0.2, n >= 30) %>%
  group_by(mu) %>%
  summarize(max_CI = max(CI_len))

# Min Interval Length For Unconstrained Analysis, n < 50
plot_df %>%
  filter(!bdd, eps == 0.2, n <= 50) %>%
  group_by(mu) %>%
  summarize(min_CI = min(CI_len))