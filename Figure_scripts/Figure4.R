source("functions.R")

# Generate Data
set.seed(6)
a <- 0; b <- 100; 
true_Y_bar <- (32.08 - a)/(b-a); true_S_sq <- (16.98/(b-a))^2; nsims <- 5e3
eps_1 <- 0.25; eps_2 <- 0.25; n <- 43; 
priv_Y_bar <- rlaplace(1, true_Y_bar, 1/(eps_1*n)); 
priv_S_sq <- rlaplace(1, true_S_sq, 1/(eps_2*n)) 

mu_0 <- (12.5 - a)/(b-a); sigma_0_sq <- (3.8/(b-a))^2
kappa_0 <- 1; nu_0 <- 1

# Run each of the Gibbs Samplers
df_bd <- Gibbs(nsims, n, priv_Y_bar, priv_S_sq, eps_1, eps_2, bounded = T, 
               flat_prior = F, mu_0 = mu_0, kappa_0 = kappa_0, nu_0 = nu_0, 
               sigma_0_sq = sigma_0_sq)
df_ub <- Gibbs(nsims, n, priv_Y_bar, priv_S_sq, eps_1, eps_2, bounded = F,
               flat_prior = F, mu_0 = mu_0, kappa_0 = kappa_0, nu_0 = nu_0, 
               sigma_0_sq = sigma_0_sq)
df_bd_flat <- Gibbs(nsims, n, priv_Y_bar, priv_S_sq, eps_1, eps_2, bounded = T, 
                    flat_prior = T)
df_ub_flat <- Gibbs(nsims, n, priv_Y_bar, priv_S_sq, eps_1, eps_2, bounded = F,
                    flat_prior = T)

# Create Shaded Region
df_fill <- data.frame(mu = b*seq(0,1,0.005), sigma_sq = 0) %>%
  mutate(ymin = 0, ymax = mu*(b-mu))

# Compute the Three Shapes
add_points_flat <- data.frame(name = "Posterior Mode", bounded = c(T, F),
                              mu = c(posterior.mode(df_bd_flat$mu, adjust = 1), 
                                     posterior.mode(df_ub_flat$mu, adjust = 1)),
                              sigma_sq = c(posterior.mode(df_bd_flat$sigma_sq, 
                                                          adjust = 1), 
                                           posterior.mode(df_ub_flat$sigma_sq, 
                                                          adjust = 1))) %>%
  rbind(data.frame(name = "Confidential Values", bounded = c(T,F), 
                   mu = true_Y_bar, sigma_sq = true_S_sq)) %>%
  rbind(data.frame(name = "Released Values", bounded = c(T,F), 
                   mu = priv_Y_bar, sigma_sq = priv_S_sq)) %>%
  mutate(bounded = factor(bounded, c(F, T), c("Unconstrained", "Constrained")),
         name = factor(name, levels = c("Released Values", "Confidential Values",
                                        "Posterior Mode", "Prior Values")),
         mu = mu*(b-a)+a, sigma_sq = sigma_sq*(b-a)^2)

# Create Figure 4
df_bd_flat %>%
  rbind(df_ub_flat) %>%
  mutate(bounded = factor(bounded, c(F, T), c("Unconstrained", "Constrained")),
         mu = mu*b, sigma_sq = sigma_sq*b^2) %>%
  ggplot(aes(x = mu, y = sigma_sq)) +
  geom_point(alpha = 0.3, size = 1) + 
  geom_ribbon(data = df_fill, aes(ymin = ymin, ymax = ymax), 
              fill = "gray", alpha = 0.5) +
  geom_point(aes(x = mu, y = sigma_sq, shape = name, color = name), 
             data = add_points_flat, size = 3) +
  facet_grid(cols = vars(bounded)) +
  scale_y_continuous(limits = c(NA, 10000)) +
  scale_shape_manual(values=c(16,17,18,15)) +
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4")) +
  labs(x = expression("Average Blood Lead ("*mu*"g/dL)"), shape = "", color = "",
       y = expression("Blood Lead Variance ("*mu*"g"^2*"/dL"^2*")")) +
  theme_tufte() +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="right",
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(-5,0,0,0),
        legend.spacing.y = unit(-0.5, "mm"))

ggsave("Figures/Figure4.pdf", width = 6, height = 2.1, dpi=600, units = "in")

## COMPUTE ADDITIONAL QUANTITIES REFERENCED IN SEC 4.2

# Proportion with Infeasible sigma_sq in Figure 3
summarize(df_ub_flat, infeasible = mean(mu > 0 & mu < 1 & sigma_sq > mu*(1-mu)))
