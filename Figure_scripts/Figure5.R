source("functions.R")

set.seed(6)
a <- 0; b <- 100; 
true_Y_bar <- (32.08 - a)/(b-a); true_S_sq <- (16.98/(b-a))^2; nsims <- 1e5
eps_1 <- 0.25; eps_2 <- 0.25; n <- 43; 
priv_Y_bar <- rlaplace(1, true_Y_bar, 1/(eps_1*n)); 
priv_S_sq <- rlaplace(1, true_S_sq, 1/(eps_2*n)) 

df_bd_flat <- Gibbs(nsims, n, priv_Y_bar, priv_S_sq, eps_1, eps_2, bounded = T, 
                    flat_prior = T)
df_ub_flat <- Gibbs(nsims, n, priv_Y_bar, priv_S_sq, eps_1, eps_2, bounded = F,
                    flat_prior = T)

df_fill_flat <- data.frame(preds = seq(a,b)) %>%
  mutate(ymin = 0, ymax = Inf)

preds_ub <- rnorm(nsims, df_ub_flat$mu, sqrt(df_ub_flat$sigma_sq))
preds_bd <- rtruncnorm(nsims, a = 0, b = 1, df_bd_flat$mu, sqrt(df_bd_flat$sigma_sq))

add_lines <- data.frame(name = "Posterior Mode", bounded = c(T, F), 
                        mu = c(posterior.mode(df_bd_flat$mu, adjust = 1), 
                               posterior.mode(df_ub_flat$mu, adjust = 1))) %>%
  rbind(data.frame(name = "Confidential Value", bounded = c(T,F), 
                   mu = true_Y_bar)) %>%
  rbind(data.frame(name = "Released Value", bounded = c(T,F), 
                   mu = priv_Y_bar)) %>%
  mutate(bounded = factor(bounded, c(F, T), c("Unconstrained", "Constrained")),
         name = factor(name, levels = c("Released Value", "Confidential Value",
                                        "Posterior Mode")),
         mu = mu*b)

data.frame(preds = c(preds_ub, preds_bd),
           bounded = c(rep("Unconstrained", nsims), 
                       rep("Constrained", nsims))) %>%
  mutate(bounded = factor(bounded, c("Unconstrained", "Constrained")), 
         preds = preds*(b-a) + a) %>%
  ggplot(aes(x = preds)) + 
  geom_histogram(binwidth = 0.1*b) +
  geom_ribbon(data = df_fill_flat, aes(ymin = ymin, ymax = ymax), 
              fill = "gray", alpha = 0.5) +
  facet_grid(~bounded) +
  scale_x_continuous(limits = c(-1*b, 2*b)) +
  labs(x = expression("New Observation's Blood Lead ("*mu*"g/dL)"),
       color = "", linetype = "", y = "Count") +
  theme_tufte() +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="right",
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(-5,0,0,0))

ggsave("Figures/Figure5.pdf", width = 6, height = 2.1, dpi=600, units = "in")


## COMPUTE ADDITIONAL QUANTITES REFERENCED IN SEC 4.2

# Proportion of negative predictions
mean(preds_ub < 0)

# Proportion of too large predictions
mean(preds_ub > 1)

# Standard deviation without constraints
sd(preds_ub)*(b-a)

# Standard deviation with constraints
sd(preds_bd)*(b-a)

# Standard deviation without constraints, but with truncation to the boundary
sd(preds_ub*(preds_ub > 0 & preds_ub < 1) + 1*(preds_ub > 1))*(b-a)

# Standard deviation without constraints, predictions from truncated Gaussian
set.seed(2)
preds_ub_trunc <- rtruncnorm(nsims, a = 0, b = 1, df_ub_flat$mu, sqrt(df_ub_flat$sigma_sq))
sd(preds_ub_trunc)*(b-a)