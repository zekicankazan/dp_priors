source("functions.R")

# Create Grid of Points
n_vec <- round(10^seq(1,3, 0.1)); eps_vec <- c(0.2, 2); 
bdd_vec <- c(FALSE, TRUE); mu_vec <- c(0.1, 0.5)
grid <- expand.grid(n = n_vec, eps = eps_vec, bdd = bdd_vec, mu = mu_vec)

# Use Saved Data to Create Dataframe For Plot
# (Saved data is created by running DCC_coverage_analysis_mode on cluster)
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