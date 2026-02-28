source("functions.R")

# Generate data
#set.seed(2)
set.seed(3)
n = 25; eps_1 <- 0.5; eps_2 <- 0.5; nsims <- 1e6
priv_Y_bar <- 1
priv_S_sq <- 0.05

mu_0 <- 0.8; sigma_0_sq <- 0.25^2
kappa_0 <- 0.5; nu_0 <- 0.5

# Run each of the Gibbs Samplers
df_bd <- Gibbs(nsims, n, priv_Y_bar, priv_S_sq, eps_1, eps_2, bounded = T, 
               flat_prior = F, mu_0 = mu_0, kappa_0 = kappa_0, nu_0 = nu_0, 
               sigma_0_sq = sigma_0_sq)
df_ub <- Gibbs(nsims, n, priv_Y_bar, priv_S_sq, eps_1, eps_2, bounded = F,
               flat_prior = F, mu_0 = mu_0, kappa_0 = kappa_0, nu_0 = nu_0, 
               sigma_0_sq = sigma_0_sq)

# Create Shaded Region
df_fill <- data.frame(mu = seq(0.5,1,0.005), sigma_sq = 0) %>%
  mutate(ymin = 0, ymax = mu*(1-mu))

df_ub_truncated <- df_ub %>%
  filter(sigma_sq <= mu*(1-mu)) %>%
  mutate(bounded = "truncated")

# Compute the Four Shapes
add_points <- data.frame(name = "Posterior Mode", bounded = c(T, "truncated", F),
                         mu = c(posterior.mode(df_bd$mu, adjust = 1), 
                                posterior.mode(df_ub_truncated$mu, adjust = 1),
                                posterior.mode(df_ub$mu, adjust = 1)),
                         sigma_sq = c(posterior.mode(df_bd$sigma_sq, adjust = 1), 
                                      posterior.mode(df_ub_truncated$sigma_sq, adjust = 1),
                                      posterior.mode(df_ub$sigma_sq, adjust = 1))) %>%
  rbind(data.frame(name = "Prior Values", bounded = c(T, "truncated", F), 
                   mu = mu_0, sigma_sq = sigma_0_sq)) %>%
  rbind(data.frame(name = "Released Values", bounded = c(T, "truncated", F), 
                   mu = priv_Y_bar, sigma_sq = priv_S_sq)) %>%
  mutate(bounded = factor(bounded, c(F, "truncated", T), 
                          c("Unconstrained", "Unconstrained & Truncated", "Constrained")),
         name = factor(name, levels = c("Released Values",
                                        "Posterior Mode", "Prior Values")),
         mu = mu, sigma_sq = sigma_sq)

add_points

plot_nsims = 1000
# # Create Figure 10
slice_sample(df_bd, n = plot_nsims) %>%
  rbind(slice_sample(df_ub_truncated, n = plot_nsims)) %>%
  rbind(slice_sample(df_ub, n = plot_nsims)) %>%
  mutate(bounded = factor(bounded, c(F, "truncated", T),
                          c("Unconstrained", "Unconstrained & Truncated", "Constrained")),
         mu = mu, sigma_sq = sigma_sq) %>%
  ggplot(aes(x = mu, y = sigma_sq)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_ribbon(data = df_fill, aes(ymin = ymin, ymax = ymax),
              fill = "gray", alpha = 0.5) +
  geom_point(aes(x = mu, y = sigma_sq, shape = name, color = name),
             data = add_points, size = 3) +
  facet_grid(cols = vars(bounded)) +
  scale_shape_manual(values=c(16,18,15)) +
  scale_color_manual(values = c("#F8766D", "#00BFC4", "#C77CFF")) +
  labs(x = expression(mu), shape = "", color = "",
       y = expression(sigma^2)) +
  theme_tufte() +
  ylim(NA, 0.4) +
  xlim(0.5, 1.3) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="right",
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(-5,0,0,0),
        legend.spacing.y = unit(-0.5, "mm"))

ggsave("Figures/Figure10.png", width = 6, height = 2.1, dpi=600, units = "in")


## COMPUTE ADDITIONAL QUANTITIES REFERENCED

# Posterior Modes
add_points %>%
  filter(name == "Posterior Mode") %>%
  mutate(sigma = sqrt(sigma_sq)) %>%
  dplyr::select(bounded, mu, sigma)

# Proportion with Infeasible sigma_sq
summarize(df_ub, infeasible = mean(mu > 0 & mu < 1 & sigma_sq > mu*(1-mu)))

# HPD Intervals for mu
round(hdi(df_bd$mu), 3)
round(hdi(df_ub_truncated$mu), 3)

# HPD Intervals for sigma_sq
round(sqrt(hdi(df_bd$sigma_sq)), 3)
round(sqrt(hdi(df_ub_truncated$sigma_sq)), 3)

