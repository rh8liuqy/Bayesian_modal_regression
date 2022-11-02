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

set.seed(2)
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

## compile stan program for normal
normal_stan_file <- "../normal_MLR.stan"
normal_stan_mod <- cmdstan_model(normal_stan_file)

## MCMC normal
normal_stan_fit <- normal_stan_mod$sample(
  data = dat,
  seed = 100,
  chains = 4,
  parallel_chains = 4,
  refresh = 0,
  iter_warmup = 5000,
  save_warmup = FALSE,
  iter_sampling = 10000,
)

## traceplot
bayesplot::mcmc_trace(normal_stan_fit$draws(c("alpha",
                                              paste0("beta[",1:P,"]"))))

## prediction interval
normal_stan_pred <- normal_stan_fit$draws(paste0("ystar[",1:N,"]"),
                                          format = "df")
normal_point_pred <- apply(normal_stan_pred, 2, median)[1:length(y)]
normal_hdi_df <- hdi(normal_stan_pred,credMass = 0.9)[,1:length(y)]
normal_df_pred <- data.frame(pred = normal_point_pred,
                             lower = normal_hdi_df[1,],
                             upper = normal_hdi_df[2,])

## model selection
normal_model_selection <- normal_stan_fit$loo()
normal_elpd_loo <- normal_model_selection$estimates[1,1]

p1 <- df1 %>%
  ggplot() +
  geom_ribbon(data = normal_df_pred, 
              aes(x = df1$X, ymin = lower, ymax=upper),
              fill = "grey70") +
  geom_line(data = normal_df_pred, 
            aes(x = df1$X, y = pred),
            color = rgb(0.8,0,0,0.8)) +
  geom_point(aes(x = X, y = y)) +
  theme_bw() +
  ylim(-60,40) +
  xlab("X") +
  ylab("y") +
  ggtitle("Likelihood: Normal (mean regression model)",
          subtitle = paste0("Coverage Rate:",
                            round(mean((y < normal_df_pred$upper) & 
                                         (y > normal_df_pred$lower)),2),
                            "||",
                            "height: ",round(mean(normal_df_pred$upper -
                                                    normal_df_pred$lower),2),
                            "||",
                            "elpd_loo: ",
                            round(normal_elpd_loo,2)))
p1

## compile stan program for median regression
median_stan_file <- "../quantile_MLR.stan"
median_stan_mod <- cmdstan_model(median_stan_file)

## MCMC median regression
median_stan_fit <- median_stan_mod$sample(
  data = dat,
  seed = 100,
  chains = 4,
  parallel_chains = 4,
  refresh = 0,
  iter_warmup = 5000,
  save_warmup = FALSE,
  iter_sampling = 10000,
)

## traceplot
bayesplot::mcmc_trace(median_stan_fit$draws(c("alpha",
                                              paste0("beta[",1:P,"]"))))

## prediction interval
median_stan_pred <- median_stan_fit$draws(paste0("ystar[",1:N,"]"),
                                          format = "df")
median_point_pred <- apply(median_stan_pred, 2, median)[1:length(y)]
median_hdi_df <- hdi(median_stan_pred,credMass = 0.9)[,1:length(y)]
median_df_pred <- data.frame(pred = median_point_pred,
                             lower = median_hdi_df[1,],
                             upper = median_hdi_df[2,])

## model selection
median_model_selection <- median_stan_fit$loo()
median_elpd_loo <- median_model_selection$estimates[1,1]

p2 <- df1 %>%
  ggplot() +
  geom_ribbon(data = median_df_pred, 
              aes(x = df1$X, ymin = lower, ymax=upper),
              fill = "grey70") +
  geom_line(data = median_df_pred, 
            aes(x = df1$X, y = pred),
            color = rgb(0.8,0,0,0.8)) +
  geom_point(aes(x = X, y = y)) +
  theme_bw() +
  ylim(-60,40) +
  xlab("X") +
  ylab("y") +
  ggtitle("Likelihood: ALD (median regression model)",
          subtitle = paste0("Coverage Rate:",
                            round(mean((y < median_df_pred$upper) & 
                                         (y > median_df_pred$lower)),2),
                            "||",
                            "height: ",round(mean(median_df_pred$upper -
                                                    median_df_pred$lower),2),
                            "||",
                            "elpd_loo: ",
                            round(median_elpd_loo,2)))
p2
grid.arrange(p1,p2,nrow=1)

## compile stan program for TPSC
TPSC_stan_file <- "../TPSC_t_MLR.stan"
TPSC_stan_mod <- cmdstan_model(TPSC_stan_file)

## MCMC median regression
TPSC_stan_fit <- TPSC_stan_mod$sample(
  data = dat,
  seed = 100,
  chains = 4,
  parallel_chains = 4,
  refresh = 0,
  iter_warmup = 5000,
  save_warmup = FALSE,
  iter_sampling = 10000,
)

## traceplot
bayesplot::mcmc_trace(TPSC_stan_fit$draws(c("alpha",
                                            paste0("beta[",1:P,"]"))))

## prediction interval
TPSC_stan_pred <- TPSC_stan_fit$draws(paste0("ystar[",1:N,"]"),
                                      format = "df")
TPSC_point_pred <- apply(TPSC_stan_pred, 2, median)[1:length(y)]
TPSC_hdi_df <- hdi(TPSC_stan_pred,credMass = 0.9)[,1:length(y)]
TPSC_df_pred <- data.frame(pred = TPSC_point_pred,
                           lower = TPSC_hdi_df[1,],
                           upper = TPSC_hdi_df[2,])

## model selection
TPSC_model_selection <- TPSC_stan_fit$loo()
TPSC_elpd_loo <- TPSC_model_selection$estimates[1,1]

p3 <- df1 %>%
  ggplot() +
  geom_ribbon(data = TPSC_df_pred, 
              aes(x = df1$X, ymin = lower, ymax=upper),
              fill = "grey70") +
  geom_line(data = TPSC_df_pred, 
            aes(x = df1$X, y = pred),
            color = rgb(0.8,0,0,0.8)) +
  geom_point(aes(x = X, y = y)) +
  theme_bw() +
  ylim(-60,40) +
  xlab("X") +
  ylab("y") +
  ggtitle("Likelihood: TPSC student t (modal regression model)",
          subtitle = paste0("Coverage Rate:",
                            round(mean((y < TPSC_df_pred$upper) & 
                                         (y > TPSC_df_pred$lower)),2),
                            "||",
                            "height: ",round(mean(TPSC_df_pred$upper -
                                                    TPSC_df_pred$lower),2),
                            "||",
                            "elpd_loo: ",
                            round(TPSC_elpd_loo,2)))
p3

pall <- grid.arrange(p1,p2,p3,ncol=1)
ggsave("left_skewed_data.pdf", pall, width = 8, height = 8)
