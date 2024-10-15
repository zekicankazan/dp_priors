library(tidyverse)
library(ggthemes)
library(purrr)

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

ggsave("Figures/Figure6.pdf", width = 6.75, height = 2, dpi=600, units = "in")