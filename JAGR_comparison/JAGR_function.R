Gibbs_JAGR <- function(nsims, n, priv_Y_bar, priv_S_sq, eps_1, eps_2, flat_prior = T, 
                       bounded = F, mu_0 = NA, kappa_0 = NA, nu_0 = NA, sigma_0_sq = NA){
  # Function to run Ju et al.'s Gibbs Sampler
  
  # Create vectors to store samples
  mu_samp <- rep(NA, nsims); sigma_sq_samp <- rep(NA, nsims)
  Y_bar_samp <- rep(NA, nsims); S_sq_samp <- rep(NA, nsims);
  
  # Initialize parameters
  Y <- runif(n); Y_bar <- mean(Y); S_sq <- var(Y)
  sigma_sq <- case_when(priv_S_sq <= 0 ~ 1e-4, priv_S_sq >= 1/4 ~ 1/4, T ~ priv_S_sq)
  accept <- rep(NA, n)
  
  
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
      sigma_sq <- 1/rgamma(1, shape = (n-2)/2, 
                           rate = ((n-1)*S_sq + n*(Y_bar - mu)^2)/2)
    }
    else if(!flat_prior & !bounded){
      sigma_sq <- 1/rgamma(1, shape = (n+nu_0)/2, 
                           rate = (nu_0*sigma_0_sq + (n-1)*S_sq + n*(Y_bar - mu)^2)/2)
    }
    else if(flat_prior & bounded){
      sigma_sq <- 1/rtrunc(1, "gamma", a = 1/(mu*(1-mu)), shape = (n-2)/2, 
                           rate = ((n-1)*S_sq + n*(Y_bar - mu)^2)/2)
    }
    else if(!flat_prior & bounded){
      sigma_sq <- 1/rtrunc(1, "gamma", a = 1/(mu*(1-mu)), shape = (n+nu_0)/2, 
                           rate = (nu_0*sigma_0_sq + (n-1)*S_sq + n*(Y_bar - mu)^2)/2)
    }
    sigma_sq_samp[t] <- sigma_sq
    
    Y_prop <- rnorm(n, mean = mu, sd = sqrt(sigma_sq))
    U <- runif(n)
    
    for(i in 1:n){
      Yi_prop <- Y_prop[i]; Yi <- Y[i]
      if(!bounded | (Yi_prop >= 0 & Yi_prop <= 1)){
        Y_bar_prop <- Y_bar + 1/n*(Yi_prop - Yi)
        S_sq_prop <- S_sq + (Yi_prop - Yi)/(n-1)*(Yi + Yi_prop - Y_bar - Y_bar_prop)
        # Formula from https://math.stackexchange.com/questions/4852929/
        
        p <- exp(-eps_1*n*(abs(priv_Y_bar - Y_bar_prop) - abs(priv_Y_bar - Y_bar))) *
          exp(-eps_2*n*(abs(priv_S_sq - S_sq_prop) - abs(priv_S_sq - S_sq)))
        if(p > U[i]){
          Y[i] <- Yi_prop; Y_bar <- Y_bar_prop; S_sq <- S_sq_prop
        }
      }
    }
    Y_bar_samp[t] <- Y_bar; S_sq_samp[t] <- S_sq
  }
  
  df <- data.frame(t = 1:nsims, mu = mu_samp, sigma_sq = sigma_sq_samp, 
                   Y_bar = Y_bar_samp, S_sq = S_sq_samp,
                   bounded = bounded, flat_prior = flat_prior)
  return(df)
}