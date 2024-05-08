#########################################
###   BASH CODE TO CALL THIS SCRIPT   ###
#########################################
#
# #!/bin/bash
# #SBATCH -e slurm.errmodule
# #SBATCH -a 1-168
# module load R/4.1.1-rhel8
# Rscript simulation_script.R
#



##############################################
###   LOAD PACKAGES AND DEFINE FUNCTIONS   ###
##############################################

suppressMessages(library(tidyverse))
suppressMessages(library(purrr))
suppressMessages(library(rmutil))
suppressMessages(library(truncdist))
suppressMessages(library(truncnorm))
suppressMessages(library(mgcv))
suppressMessages(library(HDInterval))
suppressMessages(library(MCMCglmm))

num <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## Find corresponding parameters for num
n_vec <- round(10^seq(1,3, 0.1)); eps_vec <- c(0.2, 2); 
bdd_vec <- c(FALSE, TRUE); mu_vec <- c(0.1, 0.5)
grid <- expand.grid(n = n_vec, eps = eps_vec, bdd = bdd_vec, mu = mu_vec)
n <- grid$n[num]; eps <- grid$eps[num]; bounded = grid$bdd[num]; mu = grid$mu[num]

## Set parameters
sigma_sq <- (mu*0.4)^2; eps_1 <- eps_2 <- eps/2
reps <- 1000; burn_in <- 0; nsims <- 2000

## Define necessary functions
rTGM <- function(alpha, beta, lambda, tau, upper_bound = Inf){
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



############################
###   RUN SIMULATIONS   ###
###########################


stats <- matrix(nrow = reps, ncol = 3, 
                dimnames = list(NULL, c("posterior_mode", "lower", "upper")))
set.seed(2024)
for(i in 1:reps){
  if(i %% (reps/100) == 0){print(paste0(i/reps*100, "%"))}
  stats[i,] <- run_one_experiment(i = i, mu = mu,  sigma_sq = sigma_sq, n = n, 
                                  eps_1 = eps_1,  eps_2 = eps_2, 
                                  burn_in = burn_in, nsims = nsims, 
                                  flat_prior = T, bounded = bounded)
}

## Store statistics from runs
write_csv(as.data.frame(stats),
          paste0("/coverage_analysis/n_", n, "_eps_", eps, "_bdd_", bounded, 
                 "_mu_", mu, ".csv"))
