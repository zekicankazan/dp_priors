source("functions.R")

n_vec <- round(10^seq(1,3, 0.1)); eps_vec <- c(0.2, 2); 
bdd_vec <- c(FALSE, TRUE); mu_vec <- c(0.1, 0.5)
grid <- expand.grid(n = n_vec, eps = eps_vec, bdd = bdd_vec, mu = mu_vec)
reps <- 10000; nsims <- 20000

len <- nrow(grid)
plot_df <- mutate(grid, coverage = NA, CI_len = NA, RMSE = NA)
for(num in 1:len){
  n <- grid$n[num]; eps <- grid$eps[num]; bounded <- grid$bdd[num]; mu <- grid$mu[num]
  sigma_sq <- (mu*0.4)^2; eps_1 <- eps_2 <- eps/2
  
  df <- data.frame(posterior_mode = rep(NA, reps), lower = NA, upper = NA)
  set.seed(2024)
  for(i in 1:reps){
    Y <- rtruncnorm(n, a = 0, b = 1, mean = mu, sd = sqrt(sigma_sq))
    priv_Y_bar <- rlaplace(1, mean(Y), 1/(eps_1*n))
    priv_S_sq <- rlaplace(1, var(Y), 1/(eps_2*n))
    
    out <- Gibbs(nsims, n, priv_Y_bar, priv_S_sq, eps_1, eps_2, flat_prior = T, 
                 bounded)
    
    df$posterior_mode[i] <- posterior.mode(out$mu, adjust = 1)
    df[i,c("lower", "upper")] <- hdi(out$mu, credMass = 0.95)
  }
  plot_df$CI_len[num] <- mean(df$upper - df$lower)
  plot_df$coverage[num] <- mean(mu <= df$upper & mu >= df$lower)
  plot_df$RMSE[num] <- sqrt(mean((df$posterior_mode - mu)^2))
}