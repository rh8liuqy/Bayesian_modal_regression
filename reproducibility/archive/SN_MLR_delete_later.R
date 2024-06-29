# In this R script, we explore the two parameterizations of the skew-normal distribution

rm(list = ls())
library(cmdstanr)
library(posterior)
library(bayesplot)
library(tidyverse)
library(HDInterval)

set.seed(100)
n <- 1e3

# random generation function for the univariate skew normal distribution with direct paramterization
rskewN_DP <- function(n, xi, omega, alpha) {
  U1 <- abs(rnorm(n = n))
  U2 <- rnorm(n = n)
  Z <- (alpha*U1 + U2)/sqrt(1+alpha^2)
  output <- xi + omega*Z
  return(output)
}

# random generation function for the univariate skew normal distribution with centerlized paramterization
rskewN_CP <- function(n, mu, sigma, gamma1) {
  constant_c <- (2.0 * gamma1 / (4.0 - pi))^(1.0 / 3.0)
  mu_z <- constant_c / sqrt(1.0 + constant_c^2.0)
  omega <- sigma / sqrt(1 - mu_z^2.0)
  xi <- mu - omega * mu_z
  delta <- mu_z / sqrt(2.0 / pi)
  
  if (delta > 0.0) {
    alpha <- 1.0 / sqrt(1.0 / delta^2.0 - 1.0)
  } else {
    alpha <- -1.0 / sqrt(1.0 / delta^2.0 - 1.0)
  }
  output <- rskewN_DP(n = n, xi = xi, omega = omega, alpha = alpha)
  return(output)
}

error <- rskewN_CP(n = n, mu = 0.0, sigma = 1.0, gamma1 = 0.5)
X <- matrix(runif(n,-1,1),ncol = 1)
beta <- c(-1)
eta <- 1 + as.numeric(X %*% beta)
y <- eta + error

SN_model <- cmdstan_model("./SN_MLR.stan")

SN_model_fit <- SN_model$sample(data = list(y = y, N = n, X = X, P = ncol(X)),
                                parallel_chains = 4,
                                seed = 100,
                                iter_warmup = 5000,
                                iter_sampling = 10000)
SN_model_fit$summary(c("alpha","beta"))

SN_model_fit$loo()

## prediction interval
stan_pred <- SN_model_fit$draws(paste0("ystar[",1:n,"]"),
                                format = "df")
point_pred <- apply(stan_pred, 2, median)[1:length(y)]
hdi_df <- hdi(stan_pred,credMass = 0.9)[,1:length(y)]
df_pred <- data.frame(pred = point_pred,
                      lower = hdi_df[1,],
                      upper = hdi_df[2,],
                      x = X)

df_plot <- data.frame(x = X,
                      y = y)

df_plot %>%
  ggplot() +
  geom_ribbon(data = df_pred, 
              aes(x = x, ymin = lower, ymax=upper),
              fill = "grey70") +
  geom_line(data = df_pred, 
            aes(x = x, y = pred),
            color = rgb(0.8,0,0,0.8)) +
  geom_point(aes(x = x, y = y)) +
  theme_bw()
