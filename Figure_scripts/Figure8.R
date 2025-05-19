source("simulation_study_saved.R")

# Import Saved Runs of the constrained likelihood method
grid_JAGR <- expand.grid(n = n_vec, eps = eps_vec, bdd = T, mu = mu_vec)
len <- nrow(grid_JAGR)
plot_df_JAGR <- mutate(grid_JAGR, coverage = NA, CI_len = NA, RMSE = NA)
for(num in 1:len){
  n <- grid_JAGR$n[num]; eps <- grid_JAGR$eps[num]; bdd <- T; mu <- grid_JAGR$mu[num]
  file <- paste0("JAGR_comparison/coverage_analysis_JAGR/n_", n, "_eps_", eps, "_bdd_", 
                 bdd, "_mu_", mu, "_JAGR.csv")
  if(file.exists(file)){
    df <- read.csv(file) %>%
      rename(posterior_mean = 1, lower = 2, upper = 3)
    
    plot_df_JAGR$CI_len[num] <- mean(df$upper - df$lower)
    plot_df_JAGR$coverage[num] <- mean(mu <= df$upper & mu >= df$lower)
    plot_df_JAGR$RMSE[num] <- sqrt(mean((df$posterior_mean - mu)^2))
  }
}

plot_df <- plot_df %>%
  mutate(ours = T) %>%
  rbind(mutate(plot_df_JAGR, ours = F)) %>%
  drop_na() %>%
  filter(bdd, eps == 0.2) %>%
  mutate(ours = factor(ours, c(T, F), c("Constrained Statisics", 
                                        "Constrained Likelihood")))

# Create Top Panel
p1 <- plot_df %>%
  ggplot(aes(x = n, y = CI_len, color = factor(mu), linetype = ours)) +
  geom_line(linewidth = 0.75) + geom_point(size = 1) +
  labs(color = "Truth", linetype = "",
       x = "Sample Size, n", y = "Avg. Interval Length") +
  scale_x_log10(breaks = c(10, 30, 100, 300, 1000)) +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 0.8, 0.2)) +
  scale_color_discrete(labels = c(expression("N(0.1, 0.04"^2*")"), 
                                  expression("N(0.5, 0.2"^2*")  "))) +
  scale_linetype_manual(values = c(1, 3)) +
  theme_tufte()

# Create Middle Panel
p2 <- plot_df %>%
  ggplot(aes(x = n, y = RMSE, color = factor(mu), linetype = ours)) +
  geom_line(linewidth = 0.75) + geom_point(size = 1) +
  labs(color = "Truth", linetype = "",
       x = "Sample Size, n") +
  scale_x_log10(breaks = c(10, 30, 100, 300, 1000)) +
  scale_color_discrete(labels = c(expression("N(0.1, 0.04"^2*")"), 
                                  expression("N(0.5, 0.2"^2*")  "))) +
  scale_linetype_manual(values = c(1, 3)) +
  theme_tufte() 

# Create Bottom Panel
p3 <- plot_df %>%
  ggplot(aes(x = n, y = coverage, color = factor(mu), linetype = ours)) +
  geom_line(linewidth = 0.75) + geom_point(size = 1) +
  labs(color = "Truth", linetype = "",
       x = "Sample Size, n", y = "Coverage") +
  scale_x_log10(breaks = c(10, 30, 100, 300, 1000)) +
  scale_color_discrete(labels = c(expression("N(0.1, 0.04"^2*")"), 
                                  expression("N(0.5, 0.2"^2*")  "))) +
  scale_linetype_manual(values = c(1, 3)) +
  theme_tufte()

p1/p2/p3 + plot_layout(guides = "collect", axes = "collect") & 
  theme(plot.margin=grid::unit(c(3,0,3,0), "mm"),
        panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.5))

ggsave("Figures/Figure8.pdf", width = 6, height = 4.5, dpi=600, units = "in")
