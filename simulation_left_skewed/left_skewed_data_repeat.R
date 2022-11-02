rm(list = ls())
library(cmdstanr)
library(posterior)
library(bayesplot)
library(coda)
library(MASS)
library(tidyverse)
library(lattice)
library(readxl)
library(gridExtra)
library(latex2exp)
library(parallel)
library(HDInterval)
color_scheme_set("brightblue")

simu <- function(seed,option) {
  set.seed(seed)
  ## simulate left skewed normal mixture
  N <- 30
  beta0 <- 1
  beta1 <- 1
  
  Z <- rbinom(N,1,0.05)
  epsilon <- rnorm(N,-50,1)*Z + rnorm(N,0,1)*(1-Z)
  X <- runif(N,-1,1)
  X <- as.matrix(X)
  P <- ncol(X)
  y <- beta0 + beta1*X + epsilon
  y <- as.numeric(y)
  
  df1 <- data.frame(y=y,
                    X=as.numeric(X))
  
  dat <- list(N = N,
              X = X,
              P = P,
              y = y)
  
  if (option == "normal") {
    ## compile stan program for normal
    stan_file <- "../normal_MLR.stan"
    stan_mod <- cmdstan_model(stan_file)
  }
  else if (option == "median") {
    ## compile stan program for median regression
    stan_file <- "../quantile_MLR.stan"
    stan_mod <- cmdstan_model(stan_file)
  }
  else if (option == "TPSC") {
    ## compile stan program for TPSC student t
    stan_file <- "../TPSC_t_MLR.stan"
    stan_mod <- cmdstan_model(stan_file)
  }
  ## MCMC sampling
  stan_fit <- stan_mod$sample(
    data = dat,
    seed = 100,
    chains = 4,
    parallel_chains = 4,
    refresh = 100,
    iter_warmup = 5000,
    save_warmup = FALSE,
    iter_sampling = 10000,
  )
  
  ## prediction interval
  stan_pred <- stan_fit$draws(paste0("ystar[",1:N,"]"),
                              format = "df")
  point_pred <- apply(stan_pred, 2, median)[1:length(y)]
  hdi_df <- hdi(stan_pred,credMass = 0.9)[,1:length(y)]
  df_pred <- data.frame(pred = point_pred,
                        lower = hdi_df[1,],
                        upper = hdi_df[2,])
  
  ## model selection
  model_selection <- stan_fit$loo()
  elpd_loo <- model_selection$estimates[1,1]
  df_out <- data.frame(Coverage.Rate = mean((y < df_pred$upper) & 
                                              (y > df_pred$lower)),
                       height = mean(df_pred$upper -
                                       df_pred$lower),
                       elpd_loo = elpd_loo)
  return(df_out)
}
## normal
output <- mclapply(1:300,simu,option="normal",mc.cores = 14)
saveRDS(output,"01_normal_HPD.RDS")
df_out <- bind_rows(output)
apply(df_out, 2, mean)

## ALD
output <- mclapply(1:300,simu,option="median",mc.cores = 14)
saveRDS(output,"02_quantile_HPD.RDS")
df_out <- bind_rows(output)
apply(df_out, 2, mean)

## TPSC student t
output <- mclapply(1:300,simu,option="TPSC",mc.cores = 14)
saveRDS(output,"03_TPSC_HPD.RDS")
df_out <- bind_rows(output)
apply(df_out, 2, mean)

## table creation
output <- readRDS("01_normal_HPD.RDS")
df_normal <- bind_rows(output)
output <- readRDS("02_quantile_HPD.RDS")
df_ALD <- bind_rows(output)
output <- readRDS("03_TPSC_HPD.RDS")
df_TPSC <- bind_rows(output)

df_normal$likelihood <- "normal"
df_ALD$likelihood <- "ALD"
df_TPSC$likelihood <- "TPSC student t"
library(kableExtra)
df_all <- bind_rows(df_normal,df_ALD,df_TPSC)
df_all <- df_all %>%
  group_by(likelihood) %>%
  summarise(Mean.Coverage.Rate = round(mean(Coverage.Rate),4),
            Mean.Height = round(mean(height),4),
            Mean.elpd_loo = round(mean(elpd_loo),4)
            )
kbl(df_all,
    format = "latex",
    booktabs = TRUE,
    caption = "this is caption",
    label = "this is label")
