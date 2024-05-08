library(purrr)
library(ggthemes)
library(rmutil)
library(truncdist)
library(truncnorm)
library(mgcv)
library(MASS)
library(tidyverse)
library(patchwork)
library(Matrix)
library(matrixcalc)
library(HDInterval)
library(MCMCglmm)

dTGM <- function(x, alpha, beta, lambda, tau){
  # Function to compute the PDF of the Truncated Gamma Mixture Distribution
  if(tau <= 0){
    return(dgamma(x, shape = alpha, rate = beta + lambda))
  }
  else{
    pi_1 <- exp(-lambda*tau)*pgamma(tau, alpha, beta - lambda)*gamma(alpha)/(beta - lambda)^alpha/
      (exp(-lambda*tau)*pgamma(tau, alpha, beta - lambda)*gamma(alpha)/(beta - lambda)^alpha +
         exp(lambda*tau)*pgamma(tau, alpha, beta + lambda, lower=F)*gamma(alpha)/(beta + lambda)^alpha)
    if(x <= tau){
      return(pi_1*dgamma(x, shape = alpha, rate = beta - lambda)/
               pgamma(tau, alpha, beta - lambda))
    }
    else{
      return((1-pi_1)*dgamma(x, shape = alpha, rate = beta + lambda)/
               pgamma(tau, alpha, beta + lambda, lower=F))
    }
  }
}

rTGM <- function(alpha, beta, lambda, tau, upper_bound = Inf){
  # Function to draw a random variable from a Truncated Gamma Mixture Distribution
  if(is.infinite(upper_bound) & tau <= 0){
    return(rgamma(1, shape = alpha, rate = beta + lambda))
  }
  else if(!is.infinite(upper_bound) & tau <= 0){
    return(rtrunc(1, "gamma", b = upper_bound, shape = alpha, rate = beta + lambda))
  }
  else if(tau >= upper_bound){
    return(rtrunc(1, "gamma", b = upper_bound, shape = alpha, rate = beta - lambda))
  }
  else{ #w_1, w_2 s.t. pi_1 = w_1/(w_1+w_2)
    log_w_1 <- -lambda*tau + pgamma(tau, alpha, beta - lambda, log.p = T) +
      lgamma(alpha) - alpha*log(beta - lambda)
    if(!is.infinite(upper_bound)){
      log_w_2 <- lambda*tau + log(diff(exp(pgamma(c(upper_bound, tau), alpha, 
                                                  beta + lambda, lower=F, log.p = T)))) +
        lgamma(alpha) - alpha*log(beta + lambda)
    }
    else{
      log_w_2 <- lambda*tau + pgamma(tau, alpha, beta + lambda, lower=F, log.p = T) +
        lgamma(alpha) - alpha*log(beta + lambda)
    }
    # Sample 1 with prob pi_1 = w_1/(w_1+w_2)
    ratio <- exp(log_w_1 - log_w_2) # ratio = w_1/w_2
    if(is.infinite(ratio)){ samp <- 1}
    else{ samp <- sample(c(1,2), 1, prob = c(ratio, 1))}
    if(samp == 1){
      return(rtrunc(1, "gamma", b = tau, shape = alpha, rate = beta - lambda))
    }
    else{ # samp == 2
      return(rtrunc(1, "gamma", a = tau, b = upper_bound, shape = alpha, 
                    rate = beta + lambda))
    }
  }
}

Gibbs <- function(nsims, n, priv_Y_bar, priv_S_sq, eps_1, eps_2, flat_prior = T, 
                  bounded = F, mu_0 = NA, kappa_0 = NA, nu_0 = NA, sigma_0_sq = NA){
  # Function to run our proposed Gibbs Sampler
  
  # Create vectors to store samples
  mu_samp <- rep(NA, nsims); sigma_sq_samp <- rep(NA, nsims)
  Y_bar_samp <- rep(NA, nsims); S_sq_samp <- rep(NA, nsims); 
  omega_sq_samp <- rep(NA, nsims)
  
  # Initialize parameters
  sigma_sq <- S_sq <- case_when(priv_S_sq <= 0 ~ 1e-4, priv_S_sq >= 1/4 ~ 1/4, T ~ priv_S_sq)
  Y_bar <- if_else(priv_Y_bar < 0, 0, min(priv_Y_bar, 1)); omega_sq_inv <- eps_1^2*n^2/2
  
  for(t in 1:nsims){
    # Sample mu from full conditional
    if(flat_prior & !bounded){
      mu <- rnorm(1, mean = Y_bar, sd = sqrt(sigma_sq/n))
    }
    else if(!flat_prior & !bounded){
      mu <- rnorm(1, mean = (n*Y_bar + kappa_0*mu_0)/(n+kappa_0), 
                  sd = sqrt(sigma_sq/(n+kappa_0)))
    }
    else if(sigma_sq == 1/4){
      mu <- 1/2
    }
    else if(flat_prior & bounded){
      mu <- rtruncnorm(1, a = 1/2 - sqrt(1/4 - sigma_sq), 
                       b = 1/2 + sqrt(1/4 - sigma_sq), 
                       mean = Y_bar, sd = sqrt(sigma_sq/n))
    }
    else if(!flat_prior & bounded){
      mu <- rtruncnorm(1, a = 1/2 - sqrt(1/4 - sigma_sq), 
                       b = 1/2 + sqrt(1/4 - sigma_sq), 
                       mean = (n*Y_bar + kappa_0*mu_0)/(n+kappa_0), 
                       sd = sqrt(sigma_sq/(n+kappa_0)))
    }
    mu_samp[t] <- mu
    
    # Sample sigma^2 from full conditional
    if(flat_prior & !bounded){
      sigma_sq <- 1/rtrunc(1, "gamma", a = n/(n-1)*2*eps_2, shape = (n-2)/2, 
                           rate = ((n-1)*S_sq + n*(Y_bar - mu)^2)/2)
    }
    else if(!flat_prior & !bounded){
      sigma_sq <- 1/rtrunc(1, "gamma", a = n/(n-1)*2*eps_2, shape = (n+nu_0)/2, 
                           rate = (nu_0*sigma_0_sq + (n-1)*S_sq + n*(Y_bar - mu)^2)/2)
    }
    else if(flat_prior & bounded){
      sigma_sq <- 1/rtrunc(1, "gamma", a = max(1/(mu*(1-mu)), 2*eps_2*n/(n-1)), 
                           shape = (n-2)/2, rate = ((n-1)*S_sq + n*(Y_bar - mu)^2)/2)
    }
    else if(!flat_prior & bounded){
      sigma_sq <- 1/rtrunc(1, "gamma", a = max(1/(mu*(1-mu)), 2*eps_2*n/(n-1)), 
                           shape = (n+nu_0)/2, 
                           rate = (nu_0*sigma_0_sq + (n-1)*S_sq + n*(Y_bar - mu)^2)/2)
    }
    sigma_sq_samp[t] <- sigma_sq
    
    # Sample bar{Y} from full conditional
    if(!bounded){ 
      Y_bar <- rnorm(1, sd = sqrt(1/(omega_sq_inv + n/sigma_sq)),
                     mean = (priv_Y_bar*omega_sq_inv + n*mu/sigma_sq)/
                       (omega_sq_inv + n/sigma_sq))
    }
    else{
      Y_bar <- rtruncnorm(1, a = 1/2 - sqrt(1/4 - (n-1)/n*S_sq), 
                          b = 1/2 + sqrt(1/4 - (n-1)/n*S_sq), 
                          sd = sqrt(1/(omega_sq_inv + n/sigma_sq)),
                          mean = (priv_Y_bar*omega_sq_inv + n*mu/sigma_sq)/
                            (omega_sq_inv + n/sigma_sq))
    }
    Y_bar_samp[t] <- Y_bar
    
    # Sample omega^2 from full conditional
    omega_sq_inv <- rig(1, mean = eps_1*n/abs(priv_Y_bar - Y_bar),
                        scale = 1/(eps_1^2*n^2))
    omega_sq_samp[t] <- omega_sq_inv
    
    # Sample S^2 from full conditional
    if(!bounded){
      S_sq <- rTGM(alpha = (n-1)/2, beta = (n-1)/(2*sigma_sq), lambda = eps_2*n, 
                   tau = priv_S_sq)
    }
    else{
      S_sq <- rTGM(alpha = (n-1)/2, beta = (n-1)/(2*sigma_sq), lambda = eps_2*n, 
                   tau = priv_S_sq, upper_bound = n/(n-1)*Y_bar*(1-Y_bar))
    }
    S_sq_samp[t] <- S_sq
  }
  
  df <- data.frame(t = 1:nsims, mu = mu_samp, sigma_sq = sigma_sq_samp, 
                   Y_bar = Y_bar_samp, omega_sq = omega_sq_samp, S_sq = S_sq_samp,
                   bounded = bounded, flat_prior = flat_prior)
  return(df)
}

run_one_experiment <- function(i, mu, sigma_sq, n, eps_1, eps_2, burn_in, nsims,
                               bounded = F, flat_prior = T, mu_0 = NA, 
                               kappa_0 = NA, nu_0 = NA, sigma_0_sq = NA){
  # Function to generate a single data set and run Gibbs Sampler
  
  Y <- rtruncnorm(n, a = 0, b = 1, mean = mu, sd = sqrt(sigma_sq))
  priv_Y_bar <- rlaplace(1, mean(Y), 1/(eps_1*n))
  priv_S_sq <- rlaplace(1, var(Y), 1/(eps_2*n))
  
  if(flat_prior){
    df <- Gibbs(nsims, n, priv_Y_bar, priv_S_sq, eps_1, eps_2, flat_prior = T, 
                bounded)
  }
  else{
    df <- Gibbs(nsims, n, priv_Y_bar, priv_S_sq, eps_1, eps_2, flat_prior = F, 
                bounded, mu_0, kappa_0, nu_0, sigma_0_sq)
  }
  samps <- df$mu[(burn_in+1):nsims]
  return(round(c(posterior.mode(samps, adjust = 1), hdi(samps, credMass = 0.95)),4))
}


########################################### 
###   CODE TO PRODUCE FIGURES 1 AND 3   ###
########################################### 


set.seed(6)
a <- 0; b <- 100; 
true_Y_bar <- (32.08 - a)/(b-a); true_S_sq <- (16.98/(b-a))^2; nsims <- 5e3
eps_1 <- 0.25; eps_2 <- 0.25; n <- 43; 
priv_Y_bar <- rlaplace(1, true_Y_bar, 1/(eps_1*n)); 
priv_S_sq <- rlaplace(1, true_S_sq, 1/(eps_2*n)) 

mu_0 <- (12.5 - a)/(b-a); sigma_0_sq <- (3.8/(b-a))^2
kappa_0 <- 1; nu_0 <- 1

# Run each of the four Gibbs Samplers
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

# Compute the Four Shapes for Figure 1
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

# Compute the Four Shapes for Figure 3
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
  scale_x_continuous(sec.axis = sec_axis( trans=~./b)) +
  scale_y_continuous(sec.axis = sec_axis( trans=~./b^2)) +
  scale_shape_manual(values=c(16,17,18,15)) +
  labs(x = expression("Average Blood Lead ("*mu*"g/dL)"), shape = "", color = "",
       y = expression("Blood Lead Variance ("*mu*"g"^2*"/dL"^2*")")) +
  theme_tufte() +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="bottom",
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(-5,0,0,0))

ggsave("Figures/2024_bounding_mode.png", width = 5.5, height = 3, dpi=600, 
       units = "in")

# Create Figure 3
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
  facet_grid(rows = vars(bounded)) +
  scale_x_continuous(sec.axis = sec_axis( trans=~./b)) +
  scale_y_continuous(limits = c(NA, 10000), sec.axis = sec_axis( trans=~./b^2)) +
  scale_shape_manual(values=c(16,17,18,15)) +
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4")) +
  labs(x = expression("Average Blood Lead ("*mu*"g/dL)"), shape = "", color = "",
       y = expression("Blood Lead Variance ("*mu*"g"^2*"/dL"^2*")")) +
  theme_tufte() +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="bottom",
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(-5,0,0,0))

ggsave("Figures/2024_bounding_flat_mode.png", width = 5.5, height = 3, 
       dpi=600, units = "in")

## COMPUTE ADDITIONAL QUANTITIES REFERENCED IN SEC 3.2 AND SEC 4.2

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

# Proportion with Infeasible sigma_sq in Figure 3
summarize(df_ub_flat, infeasible = mean(mu > 0 & mu < 1 & sigma_sq > mu*(1-mu)))

# HPD Intervals for mu in Figure 1
hdi(df_bd$mu)*(b-a) + a
hdi(df_ub$mu)*(b-a) + a

# HPD Intervals for sigma_sq in Figure 1
sqrt(hdi(df_bd$sigma_sq) * (b-a)^2)
sqrt(hdi(df_ub$sigma_sq) * (b-a)^2)



##################################### 
###   CODE TO PRODUCE FIGURE 4    ###
##################################### 


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
  geom_histogram(binwidth = 0.05*b) +
  geom_vline(aes(xintercept = mu, color = name, linetype = name),
             data = add_lines, linewidth = 1) +
  geom_ribbon(data = df_fill_flat, aes(ymin = ymin, ymax = ymax), 
              fill = "gray", alpha = 0.5) +
  facet_grid(bounded~ .) +
  scale_x_continuous(limits = c(-1.3*b, 2*b), 
                     sec.axis = sec_axis( trans=~./b)) +
  scale_linetype_manual(values=c("solid", "dotted", "dotdash")) +
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4")) +
  labs(x = expression("New Observation's Blood Lead Level ("*mu*"g/dL)"),
       color = "", linetype = "", y = "Count") +
  theme_tufte() +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"), legend.position="bottom",
        legend.margin=margin(0,0,0,0), legend.box.margin=margin(-5,0,0,0))

ggsave("Figures/2024_predictive.png", width = 5.5, height = 2.5, dpi=600, 
       units = "in")

## COMPUTE ADDITIONAL QUANTITES REFERENCED IN SEC 4.2

# Proportion of negative predictions
mean(preds_ub < 0)

# Proportion of too large predictions
mean(preds_ub > 1)



################################################ 
####   CODE TO PRODUCE FIGURES 2, 8 AND 9   ####
################################################

# Create Grid of Points
n_vec <- round(10^seq(1,3, 0.1)); eps_vec <- c(0.2, 2); 
bdd_vec <- c(FALSE, TRUE); mu_vec <- c(0.1, 0.5)
grid <- expand.grid(n = n_vec, eps = eps_vec, bdd = bdd_vec, mu = mu_vec)

# Use Saved Data to Create Dataframe For Plot
# (Saved data is created by running simulation_script.R on shared cluster)
len <- nrow(grid)
plot_df <- mutate(grid, coverage = NA, CI_len = NA, RMSE = NA)
for(num in 1:len){
  n <- grid$n[num]; eps <- grid$eps[num]; bdd <- grid$bdd[num]; mu <- grid$mu[num]
  file <- paste0("coverage_analysis/n_", n, "_eps_", eps, "_bdd_", 
                 bdd, "_mu_", mu, ".csv")
  if(file.exists(file)){
    df <- read.csv(file) %>%
      rename(posterior_mean = 1, lower = 2, upper = 3)
    
    plot_df$CI_len[num] <- mean(df$upper - df$lower)
    plot_df$coverage[num] <- mean(mu <= df$upper & mu >= df$lower)
    plot_df$RMSE[num] <- sqrt(mean((df$posterior_mean - mu)^2))
  }
}
plot_df <- drop_na(plot_df)

# Select Epsilon (Choose 0.2 or 2)
plot_eps <- 0.2
#plot_eps <- 2

# Create Top Panel
p1 <- plot_df %>%
  filter(n >= 30, eps == plot_eps) %>%
  mutate(bdd = factor(bdd, c(T, F))) %>%
  ggplot(aes(x = n, y = coverage, color = factor(mu), linetype = bdd)) +
  geom_line(linewidth = 0.75) + geom_point(size = 1) +
  labs(color = "Truth", linetype = "Constrained?",
       x = "Sample Size, n", y = "Coverage") +
  scale_x_log10() +
  scale_color_discrete(labels = c(expression("N(0.1, 0.04"^2*")"), 
                                  expression("N(0.5, 0.2"^2*")  "))) +
  theme_tufte() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

if(plot_eps == 2){
  p1 <- p1 + 
    geom_hline(yintercept = 0.95, linewidth = 0.75, linetype = "dotted") +
    ylim(0.89, 1)
}

# Create Middle Panel
p2 <- plot_df %>%
  filter(n >= 30, eps == plot_eps) %>%
  mutate(bdd = factor(bdd, c(T, F))) %>%
  ggplot(aes(x = n, y = CI_len, color = factor(mu), linetype = bdd)) +
  geom_line(linewidth = 0.75) + geom_point(size = 1) +
  labs(color = "Truth", linetype = "Constrained?",
       x = "Sample Size, n", y = "Interval Length") +
  scale_x_log10() +
  scale_color_discrete(labels = c(expression("N(0.1, 0.04"^2*")"), 
                                  expression("N(0.5, 0.2"^2*")  "))) +
  #scale_y_log10() +
  theme_tufte()

if(plot_eps == 0.2){
  p2 <- p2 + scale_y_continuous(breaks=seq(0, 1.5, 0.5))
}
if(plot_eps == 2){
  p2 <- p2 + scale_y_continuous(breaks=seq(0, 0.2, 0.1)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}

# Create Bottom Panel
p3 <- plot_df %>%
  filter(n >= 30, eps == plot_eps) %>%
  mutate(bdd = factor(bdd, c(T, F))) %>%
  ggplot(aes(x = n, y = RMSE, color = factor(mu), linetype = bdd)) +
  geom_line(linewidth = 0.75) + geom_point(size = 1) +
  labs(color = "Truth", linetype = "Constrained?",
       x = "Sample Size, n") +
  scale_x_log10() +
  scale_color_discrete(labels = c(expression("N(0.1, 0.04"^2*")"), 
                                  expression("N(0.5, 0.2"^2*")  "))) +
  #scale_y_log10() +
  theme_tufte()

if(plot_eps == 2){
  p3 <- p3 + scale_y_continuous(breaks=seq(0, 0.04, 0.02))
}

# Merge Panels into Single Plot
if(plot_eps == 0.2){
  p1/p2 + plot_layout(guides = "collect") & 
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
} else if(plot_eps == 2){
  p1/p2/p3 + plot_layout(guides = "collect") & 
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
}

if(plot_eps == 0.2){
  ggsave("Figures/2024_coverage_small_nolog.png", width = 5.5, height = 3.25, 
         dpi=600, units = "in")
} else if(plot_eps == 2){
  ggsave("Figures/2024_coverage_hdi_2_mode.png", width = 5.5, height = 4.5, 
         dpi=600, units = "in")
}

if(plot_eps == 0.2){
  print(p3)
  
  ggsave("Figures/2024_coverage_RMSE.png", width = 5.5, height = 2.25, 
         dpi=600, units = "in")
}

## COMPUTE ADDITIONAL QUANTITES REFERENCED IN SEC 4.1

# Linear Regression Slope of CI_length vs n
lm_CI <- plot_df %>%
  filter(eps == 0.2, !bdd) %>%
  lm(data = ., log(CI_len) ~ log(n) + mu)
lm_CI$coefficients[2] # Slope
summary(lm_CI)$r.squared # R-squared

# Linear Regression Slope of RMSE vs n
lm_RMSE <- plot_df %>%
  filter(eps == 0.2, !bdd) %>%
  lm(data = ., log(RMSE) ~ log(n) + mu)
lm_RMSE$coefficients[2] # Slope
summary(lm_RMSE)$r.squared # R-squared

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



##################################### 
###   CODE TO PRODUCE FIGURE 5    ###
##################################### 


## ADDITIONAL FUNCTIONS TO RUN SLR

Gibbs_BeSh <- function(Z, n, d, eps, Delta, eta, xi, mu_0, Lambda_0, a_0, b_0, nsims){
  # Implementation of Bernstein & Sheldon (2018)'s Gibbs-SS-Noisy Method
  
  #### Z = px1, Delta = px1, eta = dxd, xi = d*(d+1)/2 x d*(d+1)/2, mu = dx1, Lambda_0 = dxd
  #### All others constants
  p <- length(Z) # Number of releases
  
  # Create vectors to store samples
  theta_samp <- matrix(NA, nsims, d); sigma_sq_samp <- rep(NA, nsims);
  omega_sq_inv_samp <- matrix(NA, nsims, p); S_samp <- matrix(NA, nsims, p)
  
  # Initialize Parameters
  theta <- mu_0; sigma_sq <- b_0; omega_sq_inv <- (eps/p)^2/Delta^2/2
  
  for(t in 1:nsims){
    # Sample s from full Conditional
    ### Compute mu_t and Sigma_t, hardcoded for d = 2
    mu_t <- c(eta[lower.tri(eta,diag=T)], #Expectation for XTX
              theta %*% eta, #Expectation for XTY
              sigma_sq + theta %*% eta %*% theta) #Expectation for YTY
    active <- c(1:2, 4)
    Sigma_11 <- xi[active, active]
    Sigma_12 <- matrix(matrix(xi,8,2) %*% theta,4,2)[active,]
    Sigma_13 <- (xi %*% matrix(theta %*% t(theta)))[active,]
    Sigma_22 <- sigma_sq*eta + matrix(xi %*% matrix(theta %*% t(theta)),ncol = 2)
    Sigma_23 <- 2*sigma_sq*eta %*% theta  + 
      matrix(xi,2,8,T) %*% matrix(matrix(theta %*% t(theta)) %*% theta)
    Sigma_33 <- 2*sigma_sq^2 + 4*sigma_sq*theta %*% eta %*% theta +
      matrix(xi,1) %*% matrix(matrix(theta %*% t(theta)) %*% matrix(theta %*% t(theta),1))
    Sigma_t <- nearPD(rbind(cbind(Sigma_11, Sigma_12, Sigma_13),
                            cbind(t(Sigma_12), Sigma_22, Sigma_23),
                            c(Sigma_13, t(Sigma_23), Sigma_33)))$mat
    
    ### Sample S
    Sigma_3 <- solve(solve(Sigma_t)/n + diag(omega_sq_inv))
    mu_3 <- Sigma_3 %*% (solve(Sigma_t) %*% mu_t + diag(omega_sq_inv) %*% Z)
    S <- mvrnorm(1, mu = mu_3, Sigma = Sigma_3)
    S_samp[t,] <- as.vector(S)
    
    # Sample theta, sigma_sq from full Conditional
    ### Obtain current sufficient statistics
    XTX <- matrix(NA,d,d); XTX[lower.tri(t(XTX), diag = T)] <- S[1:(d*(d+1)/2)]
    XTX[upper.tri(XTX)] <- (t(XTX))[upper.tri(XTX)]
    XTY <- matrix(S[(d*(d+1)/2 + 1):(d*(d+1)/2 + d)],d,1); YTY <- S[length(S)]
    
    B <- as.matrix(nearPD(rbind(cbind(XTX, XTY), cbind(t(XTY), YTY)))$mat)
    XTX <- B[1:d, 1:d]; XTY <- matrix(B[1:d, d+1]); YTY <- B[d+1, d+1]
    
    ### Compute hyperparameters
    mu_n <- solve(XTX + Lambda_0) %*% (XTY + Lambda_0 %*% mu_0)
    Lambda_n <- XTX + Lambda_0
    a_n <- a_0 + n/2
    b_n <- b_0 + (YTY + t(mu_0) %*% Lambda_0 %*% mu_0 - 
                    as.numeric(t(mu_n) %*% Lambda_n %*% mu_n))/2
    
    ### Sample
    sigma_sq <- 1/rgamma(1, shape = a_n,  rate = b_n)
    theta <- mvrnorm(1, mu = mu_n, Sigma = solve(nearPD(Lambda_n)$mat)*sigma_sq)
    
    ### Save samples
    sigma_sq_samp[t] <- sigma_sq
    theta_samp[t,] <- theta
    
    # Sample omega_sq_inv from full conditional
    omega_sq_inv <- rig(p, mean = as.vector(eps/p/Delta/abs(Z - S)), 
                        scale = Delta^2/(eps/p)^2)
    omega_sq_inv_samp[t,] <- omega_sq_inv
  }
  return(list(t = 1:t, theta = theta_samp, sigma_sq = sigma_sq_samp, 
              omega_sq_inv = omega_sq_inv_samp, S = S_samp))
}

Gibbs_BeSh_bound <- function(Z, n, d, eps, Delta, eta, xi, mu_0, Lambda_0, a_0, b_0, nsims){
  # Implementation of Gibbs-SS-Noisy with Constraints
  
  #### Z = px1, Delta = px1, eta = dxd, xi = d*(d+1)/2 x d*(d+1)/2, mu = dx1, Lambda_0 = dxd
  #### All others constants
  p <- length(Z) # Number of releases
  
  # Create vectors to store samples
  theta_samp <- matrix(NA, nsims, d); sigma_sq_samp <- rep(NA, nsims);
  omega_sq_inv_samp <- matrix(NA, nsims, p); S_samp <- matrix(NA, nsims, p)
  
  # Initialize Parameters
  theta <- mu_0; sigma_sq <- b_0; omega_sq_inv <- (eps/p)^2/Delta^2/2
  
  for(t in 1:nsims){
    #if(t %% 10 == 0){print(paste0("t = ", t, "  ", Sys.time()))}
    
    # Sample s from full Conditional
    ### Compute mu_t and Sigma_t, hardcoded for d = 2
    mu_t <- c(eta[,2], #Expectation for XTX
              theta %*% eta, #Expectation for XTY
              sigma_sq + theta %*% eta %*% theta) #Expectation for YTY
    active <- 3:4
    Sigma_11 <- xi[active, active]
    Sigma_12 <- matrix(matrix(xi,8,2) %*% theta,4,2)[active,]
    Sigma_13 <- (xi %*% matrix(theta %*% t(theta)))[active,]
    Sigma_22 <- sigma_sq*eta + matrix(xi %*% matrix(theta %*% t(theta)),ncol = 2)
    Sigma_23 <- 2*sigma_sq*eta %*% theta  + 
      matrix(xi,2,8,T) %*% matrix(matrix(theta %*% t(theta)) %*% theta)
    Sigma_33 <- 2*sigma_sq^2 + 4*sigma_sq*theta %*% eta %*% theta +
      matrix(xi,1) %*% matrix(matrix(theta %*% t(theta)) %*% matrix(theta %*% t(theta),1))
    Sigma_t <- nearPD(rbind(cbind(Sigma_11, Sigma_12, Sigma_13),
                            cbind(t(Sigma_12), Sigma_22, Sigma_23),
                            c(Sigma_13, t(Sigma_23), Sigma_33)))$mat
    
    ### Sample S
    Sigma_3 <- solve(solve(Sigma_t)/n + diag(omega_sq_inv))
    mu_3 <- Sigma_3 %*% (solve(Sigma_t) %*% mu_t + diag(omega_sq_inv) %*% Z)
    
    flag <- F
    try <- 0
    while(!flag){
      try <- try + 1
      S <- mvrnorm(1, mu = mu_3, Sigma = Sigma_3)
      
      B <- matrix(c(n, S[c(1,3,1,2,4,3:5)]),3)
      flag <- (S[1] >= 0 & S[1] <= n) & (S[2] >= 0 & S[2] <= S[1]) &
        (S[3] >= 0 & S[3] <= n) & (S[4] >= 0 & S[4] <= S[1] & S[4] <= S[3]) &
        (S[5] >= 0 & S[5] <= S[3]) & is.positive.definite(B)
    }
    S_samp[t,] <- as.vector(S)
    
    # Sample theta, sigma_sq from full Conditional
    ### Obtain current sufficient statistics
    XTX <- matrix(NA,d,d); XTX[lower.tri(t(XTX), diag = T)] <- c(n, S[1:(d*(d+1)/2 - 1)])
    XTX[upper.tri(XTX)] <- (t(XTX))[upper.tri(XTX)]
    XTY <- matrix(S[(d*(d+1)/2):(d*(d+1)/2 + d - 1)],d,1); YTY <- S[length(S)]
    
    #B <- as.matrix(nearPD(rbind(cbind(XTX, XTY), cbind(t(XTY), YTY)))$mat)
    #XTX <- B[1:d, 1:d]; XTY <- matrix(B[1:d, d+1]); YTY <- B[d+1, d+1]
    
    ### Compute hyperparameters
    mu_n <- solve(XTX + Lambda_0) %*% (XTY + Lambda_0 %*% mu_0)
    Lambda_n <- XTX + Lambda_0
    a_n <- a_0 + n/2
    b_n <- b_0 + (YTY + t(mu_0) %*% Lambda_0 %*% mu_0 - 
                    as.numeric(t(mu_n) %*% Lambda_n %*% mu_n))/2
    
    ### Sample
    flag = F
    while(!flag){
      sigma_sq <- 1/rtrunc(1, "gamma", b = 1/4, shape = a_n,  rate = b_n)
      theta <- mvrnorm(1, mu = mu_n, Sigma = solve(nearPD(Lambda_n)$mat)*sigma_sq)
      flag <- (theta[1] >= 0 & theta[1] <= 1) | (theta[1] > 1 & theta[2] <= 1 - theta[1]) |
        (theta[1] < 0 & theta[2] >= -theta[1])
    }
    
    ### Save samples
    sigma_sq_samp[t] <- sigma_sq
    theta_samp[t,] <- theta
    
    # Sample omega_sq_inv from full conditional
    omega_sq_inv <- rig(p, mean = as.vector(eps/p/Delta/abs(Z - S)), 
                        scale = Delta^2/(eps/p)^2)
    omega_sq_inv_samp[t,] <- omega_sq_inv
  }
  return(list(t = 1:t, theta = theta_samp, sigma_sq = sigma_sq_samp, 
              omega_sq_inv = omega_sq_inv_samp, S = S_samp))
}

# Import Data
x20 <- read.csv("https://people.sc.fsu.edu/~jburkardt/datasets/regression/x20.txt", 
                sep = "", header = F, comment.char = "#",
                col.names = c("ID", "Intercept", "pop_size", "births", "wine", "liquor",
                              "death")) %>%
  tail(-9) %>%
  remove_rownames() %>%
  column_to_rownames(var = "ID") %>%
  mutate(across(.fns = as.numeric))

# Compute X and Y
n <- nrow(x20)
Y <- (x20$death - min(x20$death))/(max(x20$death) - min(x20$death))
drink <- x20$wine + x20$liquor; x1 <- (drink - min(drink))/(max(drink) - min(drink))
X <- matrix(c(x20$Intercept, x1), ncol = 2)

# Compute Noisy Statistics
set.seed(10)
eps <- 0.1*11; n <- nrow(X); d <- ncol(X); p <- d*(d+1)/2 + d + 1
Delta <- rep(1, p); nsims <- 5000
mu_0 <- c(1,0); Lambda_0 <- diag(c(.25, .25)); a_0 <- 20; b_0 <- 0.5

XTX_true <- t(X) %*% X; XTY_true <- t(X) %*% Y; YTY_true <- t(Y) %*% Y
S_true <- c(XTX_true[lower.tri(XTX_true,diag=T)], XTY_true, YTY_true)
Z <- rlaplace(p, S_true, Delta/(eps/p))

X_4th_priv <- rlaplace(n = 5, m = map_dbl(.x = 0:4, .f = function(i){sum(X[,2]^i)}),
                       s = 1/(eps/11))

# Compute Quantities from Bernstein & Sheldon (2018)
eta <- matrix(X_4th_priv[c(1,2,2,3)],2,2)/n
eta_ijkl <- matrix(X_4th_priv[c(1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5)], 4, 4)/n
xi <- eta_ijkl - as.vector(eta) %*% t(as.vector(eta))

# Run Both Algorithms
l = Gibbs_BeSh(Z, n, d, eps*6/11, Delta, eta, xi, mu_0, Lambda_0, a_0, b_0, nsims)
l_bd = Gibbs_BeSh_bound(Z[2:length(Z)], n, d, eps*6/11, Delta[2:length(Z)], eta, xi, 
                        mu_0, Lambda_0, a_0, b_0, nsims) 

# Create Shaded Region for Left Panels
df_fill_dat <- expand.grid(YT1 = seq(0, 42, 1), YTY = seq(0, 42, 1)) %>%
  mutate(ymax = YT1*(YT1 >= 0 & YT1 <= n), ymin = 0)

# Find Infeasible Points in Left Panels
df_infeasible_dat <- as.data.frame(l$S[,2:6]) %>%
  mutate(bdd = F) %>% 
  rename(X1T1 = V1, X1TX1 = V2, YT1 = V3, X1TY = V4, YTY = V5) %>%
  mutate(bdd = factor(bdd, c(F, T), c("Unconstrained", "Constrained"))) %>%
  filter(!(0 < YTY & YTY < YT1 & YT1 < n))

# Create Left Panels
p1 <- rbind(mutate(as.data.frame(l$S[,2:6]), bdd = F),
            mutate(as.data.frame(l_bd$S), bdd = T)) %>%
  rename(X1T1 = V1, X1TX1 = V2, YT1 = V3, X1TY = V4, YTY = V5) %>%
  mutate(bdd = factor(bdd, c(F, T), c("Unconstrained", "Constrained"))) %>%
  ggplot(aes(x = YT1, y = YTY)) + 
  geom_point(size = 1, alpha = 0.3) +
  geom_point(data = df_infeasible_dat, aes(x = YT1, y = YTY), color = "#F8766D",
             size = 1, alpha = 0.5) +
  facet_grid(bdd~.) +
  geom_ribbon(data = df_fill_dat, aes(ymin = ymin, ymax = ymax), fill = "gray", 
              alpha = 0.5) +
  labs(x = expression("Y"^"T"*"1"), y = expression("Y"^"T"*"Y")) +
  theme_tufte() +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank())

# Create Shaded Region for Right Panels
df_fill_pars <- expand.grid(Slope = seq(-4, 4, 0.04), Intercept = seq(-5,5, 0.05)) %>%
  mutate(ymax = 1*(Slope >= 0) + (1-Slope)*(Slope < 0),
         ymin = (-Slope)*(Slope >= 0))

# Find Infeasible Points in Right Panels
df_infeasible_pars <- as.data.frame(l$theta) %>%
  mutate(bdd = F) %>% 
  rename(Intercept = V1, Slope = V2) %>%
  mutate(bdd = factor(bdd, c(F, T), c("Unconstrained", "Constrained"))) %>%
  filter((Slope >= 0) & !(Intercept > -Slope & Intercept < 1) |
           (Slope < 0) & !(Intercept > 0 & Intercept < 1-Slope))

# Create Right Panels
p2 <- rbind(mutate(as.data.frame(l$theta), bdd = F), 
            mutate(as.data.frame(l_bd$theta), bdd = T)) %>%
  rename(Intercept = V1, Slope = V2) %>%
  mutate(bdd = factor(bdd, c(F, T), c("Unconstrained", "Constrained"))) %>%
  ggplot(aes(x = Slope, y = Intercept)) +
  geom_point(size = 1, alpha = 0.3) +
  geom_point(data = df_infeasible_pars,aes(x = Slope, y = Intercept), color = "#F8766D",
             size = 1, alpha = 0.5) +
  facet_grid(bdd~.) +
  xlim(-4,4) +
  labs(x = expression(theta[1]), y = expression(theta[0])) +
  geom_ribbon(data = df_fill_pars, aes(ymin = ymin, ymax = ymax), fill = "gray", alpha = 0.5) + 
  theme_tufte()

# Create Figure 5
p1 + p2

ggsave("Figures/2024_BeSh_comparison.png", width = 5.5, height = 3, dpi=600, units = "in")

## COMPUTE ADDITIONAL QUANTITES REFERENCED IN SEC 5

# Proportion Infeasible for Left Panel
nrow(df_infeasible_dat)/nsims

# Proportion Infeasible for Right Panel
nrow(df_infeasible_pars)/nsims



##################################### 
###   CODE TO PRODUCE FIGURE 6    ###
##################################### 


data.frame(x = seq(0.01, 4, 0.01)) %>%
  mutate(Gamma = dgamma(x, shape = 2, rate = 2),
         TGM = map_dbl(.x = x, .f = dTGM, alpha = 2, beta = 2, lambda = 1, tau = 1)) %>%
  pivot_longer(-x, names_to = "Distribution", values_to = "Density") %>%
  mutate(Distribution = factor(Distribution, levels = c("TGM", "Gamma"))) %>%
  ggplot(aes(x = x, y = Density, linetype = Distribution)) + 
  labs(x = "", linetype = "") +
  geom_line(linewidth = 1) +
  theme_tufte() +
  theme(legend.position=c(.85,.5)) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

ggsave("Figures/2024_TGM.png", width = 5.5, height = 2, dpi=600, units = "in")