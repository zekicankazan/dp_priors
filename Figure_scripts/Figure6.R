source("functions.R")
source("JAGR_comparison/JAGR_function.R")

set.seed(5)
true_Y_bar <- .5; true_S_sq <- 0.2^2; nsims <- 2.5e4
eps_1 <- 0.25; eps_2 <- 0.25; n <- 50; 
priv_Y_bar <- rlaplace(1, true_Y_bar, 1/(eps_1*n)); 
priv_S_sq <- rlaplace(1, true_S_sq, 1/(eps_2*n)) 

# Run each sampler
df_ours <- Gibbs(nsims, n, priv_Y_bar, priv_S_sq, eps_1, eps_2, 
                 bounded = T, flat_prior = T)

df_JAGR <- Gibbs_JAGR(nsims, n, priv_Y_bar, priv_S_sq, eps_1, eps_2, 
                      bounded = T, flat_prior = T)

# Create feasible region
df_fill <- data.frame(mu = seq(0,1,0.005), sigma_sq = 0) %>%
  mutate(ymin = 0, ymax = mu*(1-mu))

# Compute additional points
add_points_comp <- data.frame(name = "Posterior Mode", ours = c(T, F),
                              mu = c(posterior.mode(df_ours$mu, adjust = 1), 
                                     posterior.mode(df_JAGR$mu, adjust = 1)),
                              sigma_sq = c(posterior.mode(df_ours$sigma_sq, 
                                                          adjust = 1), 
                                           posterior.mode(df_JAGR$sigma_sq, 
                                                          adjust = 1))) %>%
  rbind(data.frame(name = "Confidential Values", ours = c(T,F), 
                   mu = true_Y_bar, sigma_sq = true_S_sq)) %>%
  rbind(data.frame(name = "Released Values", ours = c(T,F), 
                   mu = priv_Y_bar, sigma_sq = priv_S_sq)) %>%
  mutate(ours = factor(ours, c(T, F), c("Constrained Statisics", "Constrained Likelihood")),
         name = factor(name, levels = c("Released Values", "Confidential Values",
                                        "Posterior Mode")))

# Create Figure 6
mutate(df_ours, ours = T) %>%
  dplyr::select(-omega_sq) %>%
  rbind(mutate(df_JAGR, ours = F)) %>%
  filter(t %% 5 == 0) %>%
  mutate(ours = factor(ours, c(T, F), c("Constrained Statisics", "Constrained Likelihood"))) %>%
  ggplot(aes(x = mu, y = sigma_sq)) +
  geom_point(alpha = 0.3, size = 1) + 
  geom_ribbon(data = df_fill, aes(ymin = ymin, ymax = ymax), 
              fill = "gray", alpha = 0.5) +
  geom_point(aes(x = mu, y = sigma_sq, shape = name, color = name), 
             data = add_points_comp, size = 3) +
  facet_grid(cols = vars(ours)) +
  scale_shape_manual(values=c(16,17,18,15)) +
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4")) +
  labs(x = expression(mu), shape = "", color = "", y = expression(sigma^2)) +
  theme_tufte() +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

ggsave("Figures/Figure6.pdf", width = 6, height = 2.1, dpi=600, units = "in")
