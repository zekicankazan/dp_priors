source("functions.R")

# Generate data
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

# Create Shaded Region
df_fill <- data.frame(mu = b*seq(0,1,0.005), sigma_sq = 0) %>%
  mutate(ymin = 0, ymax = mu*(b-mu))

# Compute the Four Shapes
add_points <- data.frame(name = "Posterior Mode", bounded = c(T, F),
                         mu = c(posterior.mode(df_bd$mu, adjust = 1), 
                                posterior.mode(df_ub$mu, adjust = 1)),
                         sigma_sq = c(posterior.mode(df_bd$sigma_sq, adjust = 1), 
                                      posterior.mode(df_ub$sigma_sq, adjust = 1))) %>%
  rbind(data.frame(name = "Prior Values", bounded = c(T,F), 
                   mu = mu_0, sigma_sq = sigma_0_sq)) %>%
  rbind(data.frame(name = "Confidential Values", bounded = c(T,F), 
                   mu = true_Y_bar, sigma_sq = true_S_sq)) %>%
  rbind(data.frame(name = "Released Values", bounded = c(T,F), 
                   mu = priv_Y_bar, sigma_sq = priv_S_sq)) %>%
  mutate(bounded = factor(bounded, c(F, T), c("Unconstrained", "Constrained")),
         name = factor(name, levels = c("Released Values", "Confidential Values",
                                        "Posterior Mode", "Prior Values")),
         mu = mu*b, sigma_sq = sigma_sq*b^2)

# Create Figure 1
df_bd %>%
  rbind(df_ub) %>%
  mutate(bounded = factor(bounded, c(F, T), c("Unconstrained", "Constrained")),
         mu = mu*b, sigma_sq = sigma_sq*b^2) %>%
  ggplot(aes(x = mu, y = sigma_sq)) +
  geom_point(alpha = 0.3, size = 1) + 
  geom_ribbon(data = df_fill, aes(ymin = ymin, ymax = ymax), 
              fill = "gray", alpha = 0.5) +
  geom_point(aes(x = mu, y = sigma_sq, shape = name, color = name),
             data = add_points, size = 3) +
  facet_grid(rows = vars(bounded)) +
  scale_shape_manual(values=c(16,17,18,15)) +
  labs(x = expression("Average Blood Lead ("*mu*"g/dL)"), shape = "", color = "",
       y = expression("Blood Lead Variance ("*mu*"g"^2*"/dL"^2*")")) +
  theme_tufte() +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="bottom",
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(-5,0,0,0),
        legend.spacing.y = unit(-0.5, "mm")) +
  guides(colour=guide_legend(ncol=2,nrow=2,byrow=TRUE),
         shape = guide_legend(ncol=2,nrow=2,byrow=TRUE))

ggsave("Figures/Figure1.pdf", width = 3, height = 3.4, dpi=600, units = "in")


## COMPUTE ADDITIONAL QUANTITIES REFERENCED IN SEC 3.2

# Released Statistics
priv_Y_bar*(b-a) + a
sqrt(priv_S_sq)*(b-a)

# Posterior Modes
add_points %>%
  filter(name == "Posterior Mode") %>%
  mutate(sigma = sqrt(sigma_sq)) %>%
  dplyr::select(bounded, mu, sigma)

# Proportion with Infeasible sigma_sq in Figure 1
summarize(df_ub, infeasible = mean(mu > 0 & mu < 1 & sigma_sq > mu*(1-mu)))

# HPD Intervals for mu in Figure 1
hdi(df_bd$mu)*(b-a) + a
hdi(df_ub$mu)*(b-a) + a

# HPD Intervals for sigma_sq in Figure 1
sqrt(hdi(df_bd$sigma_sq) * (b-a)^2)
sqrt(hdi(df_ub$sigma_sq) * (b-a)^2)