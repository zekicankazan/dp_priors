source("simulation_study_saved.R")

plot_eps <- 2

# Create Top Panel
p1 <- plot_df %>%
  filter(n >= 30, eps == plot_eps) %>%
  mutate(bdd = factor(bdd, c(T, F))) %>%
  ggplot(aes(x = n, y = CI_len, color = factor(mu), linetype = bdd)) +
  geom_line(linewidth = 0.75) + geom_point(size = 1) +
  labs(color = "Truth", linetype = "Constrained?", x = "Sample Size, n", y = "Avg. Interval Length") +
  scale_x_log10() +
  scale_y_continuous(limits = c(NA, 0.3), breaks=seq(0, 0.3, 0.1)) +
  scale_color_discrete(labels = c(expression("N(0.1, 0.04"^2*")"), expression("N(0.5, 0.2"^2*")  "))) +
  theme_tufte()  +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

# Create Middle Panel
p2 <- plot_df %>%
  filter(n >= 30, eps == plot_eps) %>%
  mutate(bdd = factor(bdd, c(T, F))) %>%
  ggplot(aes(x = n, y = RMSE, color = factor(mu), linetype = bdd)) +
  geom_line(linewidth = 0.75) + geom_point(size = 1) +
  labs(color = "Truth", linetype = "Constrained?",
       x = "Sample Size, n") +
  scale_x_log10() +
  scale_y_continuous(limits = c(NA, 0.06), breaks=seq(0, 0.06, 0.02)) +
  scale_color_discrete(labels = c(expression("N(0.1, 0.04"^2*")"), 
                                  expression("N(0.5, 0.2"^2*")  "))) +
  theme_tufte() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Create Bottom Panel
p3 <- plot_df %>%
  filter(n >= 30, eps == plot_eps) %>%
  mutate(bdd = factor(bdd, c(T, F))) %>%
  ggplot(aes(x = n, y = coverage, color = factor(mu), linetype = bdd)) +
  geom_line(linewidth = 0.75) + geom_point(size = 1) +
  labs(color = "Truth", linetype = "Constrained?", x = "Sample Size, n", y = "Coverage") +
  scale_x_log10() +
  scale_y_continuous(limits = c(0.9, 1), breaks=seq(0.9, 1, 0.025)) +
  scale_color_discrete(labels = c(expression("N(0.1, 0.04"^2*")"), expression("N(0.5, 0.2"^2*")  "))) +
  theme_tufte()

# Create Figure 9
(p1/p2/p3) + plot_layout(guides = "collect") & 
  theme(plot.margin=grid::unit(c(3,0,3,0), "mm"),
        panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.5))

ggsave("Figures/Figure9.pdf", width = 6.75, height = 4.5, dpi=600, units = "in")
