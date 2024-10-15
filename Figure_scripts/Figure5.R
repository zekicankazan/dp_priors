library(tidyverse)
library(purrr)
library(ggthemes)
library(rmutil)
library(truncdist)
library(patchwork)
library(Matrix)
library(matrixcalc)
library(MASS)
library(mgcv)

## ADDITIONAL FUNCTION TO RUN SLR
Gibbs_BeSh <- function(Z, n, d, eps, Delta, eta, xi, mu_0, Lambda_0, a_0, b_0, nsims,
                       bounded = F){
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
    if(!bounded){
      mu_t <- c(eta[lower.tri(eta,diag=T)], #Expectation for XTX
                theta %*% eta, #Expectation for XTY
                sigma_sq + theta %*% eta %*% theta) #Expectation for YTY
      active <- c(1:2, 4)
    }
    else if(bounded){
      mu_t <- c(eta[,2], #Expectation for XTX
                theta %*% eta, #Expectation for XTY
                sigma_sq + theta %*% eta %*% theta) #Expectation for YTY
      active <- 3:4
    }
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
    if(!bounded){
      S <- mvrnorm(1, mu = mu_3, Sigma = Sigma_3)
    }
    else if(bounded){
      flag <- F
      while(!flag){
        S <- mvrnorm(1, mu = mu_3, Sigma = Sigma_3)
        B <- matrix(c(n, S[c(1,3,1,2,4,3:5)]),3)
        flag <- (S[1] >= 0 & S[1] <= n) & (S[2] >= 0 & S[2] <= S[1]) &
          (S[3] >= 0 & S[3] <= n) & (S[4] >= 0 & S[4] <= S[1] & S[4] <= S[3]) &
          (S[5] >= 0 & S[5] <= S[3]) & is.positive.definite(B)
      }}
    S_samp[t,] <- as.vector(S)
    
    # Sample theta, sigma_sq from full Conditional
    ### Obtain current sufficient statistics
    XTX <- matrix(NA,d,d); 
    if(!bounded){
      XTX[lower.tri(t(XTX), diag = T)] <- S[1:(d*(d+1)/2)]
      XTY <- matrix(S[(d*(d+1)/2 + 1):(d*(d+1)/2 + d)],d,1)
    }
    else if(bounded){
      XTX[lower.tri(t(XTX), diag = T)] <- c(n, S[1:(d*(d+1)/2 - 1)])
      XTY <- matrix(S[(d*(d+1)/2):(d*(d+1)/2 + d - 1)],d,1)
    }
    XTX[upper.tri(XTX)] <- (t(XTX))[upper.tri(XTX)]
    YTY <- S[length(S)]
    
    if(!bounded){
      B <- as.matrix(nearPD(rbind(cbind(XTX, XTY), cbind(t(XTY), YTY)))$mat)
      XTX <- B[1:d, 1:d]; XTY <- matrix(B[1:d, d+1]); YTY <- B[d+1, d+1]
    }
    
    ### Compute hyperparameters
    mu_n <- solve(XTX + Lambda_0) %*% (XTY + Lambda_0 %*% mu_0)
    Lambda_n <- XTX + Lambda_0
    a_n <- a_0 + n/2
    b_n <- b_0 + (YTY + t(mu_0) %*% Lambda_0 %*% mu_0 - 
                    as.numeric(t(mu_n) %*% Lambda_n %*% mu_n))/2
    
    ### Sample
    if(!bounded){
      sigma_sq <- 1/rgamma(1, shape = a_n,  rate = b_n)
      theta <- mvrnorm(1, mu = mu_n, Sigma = solve(nearPD(Lambda_n)$mat)*sigma_sq)
    }
    else if(bounded){
      flag = F
      while(!flag){
        sigma_sq <- 1/rtrunc(1, "gamma", b = 1/4, shape = a_n,  rate = b_n)
        theta <- mvrnorm(1, mu = mu_n, Sigma = solve(nearPD(Lambda_n)$mat)*sigma_sq)
        flag <- (theta[1] >= 0 & theta[1] <= 1) | (theta[1] > 1 & theta[2] <= 1 - theta[1]) |
          (theta[1] < 0 & theta[2] >= -theta[1])
      }}
    
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
#l_bd = Gibbs_BeSh_bound(Z[2:length(Z)], n, d, eps*6/11, Delta[2:length(Z)], eta, xi, 
#                        mu_0, Lambda_0, a_0, b_0, nsims) 
l_bd = Gibbs_BeSh(Z[2:length(Z)], n, d, eps*6/11, Delta[2:length(Z)], eta, xi, 
                  mu_0, Lambda_0, a_0, b_0, nsims, bounded = T)

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
  scale_y_continuous(breaks = seq(0, 90, 30)) +
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
p1 + p2 & 
  theme(plot.margin=grid::unit(c(0.5,0,0,0), "mm"))

ggsave("Figures/Figure5.pdf", width = 6.75, height = 2.15, dpi=600, units = "in")



## COMPUTE ADDITIONAL QUANTITES REFERENCED IN SEC 5

# Proportion Infeasible for Left Panel
nrow(df_infeasible_dat)/nsims

# Proportion Infeasible for Right Panel
nrow(df_infeasible_pars)/nsims