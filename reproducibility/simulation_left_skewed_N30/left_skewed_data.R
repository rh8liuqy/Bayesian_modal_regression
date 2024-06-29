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
library(kableExtra)
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

par.est <- normal_stan_fit$summary(c("alpha",paste0("beta[",1:P,"]")))
print(par.est)

## traceplot
pdf("traceplot_normal.pdf",height = 3,width = 8)
bayesplot::mcmc_trace(normal_stan_fit$draws(c("alpha",
                                              paste0("beta[",1:P,"]"),
                                              "sigma")),
                      facet_args = list(ncol = 2))
dev.off()

## summary statistics
normal_stan_fit$summary(c("alpha",
                          paste0("beta[",1:P,"]"),
                          "sigma")) %>%
  kbl(format = "latex", booktabs = TRUE, digits = 2) %>%
  kable_styling()

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
  ggtitle("Likelihood: Normal (classic mean regression model)",
          subtitle = paste0("Coverage Rate:",
                            round(mean((y < normal_df_pred$upper) & 
                                         (y > normal_df_pred$lower)),2),
                            "||",
                            "Width: ",round(mean(normal_df_pred$upper -
                                                    normal_df_pred$lower),2),
                            "||",
                            "ELPD: ",
                            round(normal_elpd_loo,2)))
p1

## true error density
df_error <- data.frame(x = seq(from = -60, to = 40, length.out = 1000))
true_error_density <- function(x) {
  p1 <- 0.05*dnorm(x = x, mean = -50, sd = 1)
  p2 <- 0.95*dnorm(x = x, mean = 0, sd = 1)
  output <- p1 + p2
  return(output)
}
df_error$density <- true_error_density(x = df_error$x)
df_error$loglik <- "True"

## calculation of residuals
sigma_summary <- normal_stan_fit$summary("sigma")
df_normal_residual <- data.frame(x = seq(from = -60, to = 40, length.out = 1000))
df_normal_residual$density <- dnorm(x = df_normal_residual$x,
                                    mean = 0,
                                    sd = sigma_summary$mean)
df_normal_residual$loglik <- "Normal"

df_normal_residual %>%
  ggplot(aes(x = x, y = density)) +
  geom_line()

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

par.est <- median_stan_fit$summary(c("alpha",
                                     paste0("beta[",1:P,"]"),
                                     "sigma"))
print(par.est)

## traceplot
pdf("traceplot_ALD.pdf",height = 3,width = 8)
bayesplot::mcmc_trace(median_stan_fit$draws(c("alpha",
                                              paste0("beta[",1:P,"]"),
                                              "sigma")),
                      facet_args = list(ncol = 2))
dev.off()

## summary statistics
median_stan_fit$summary(c("alpha",
                          paste0("beta[",1:P,"]"),
                          "sigma")) %>%
  kbl(format = "latex", booktabs = TRUE, digits = 2) %>%
  kable_styling()

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
                            "Width: ",round(mean(median_df_pred$upper -
                                                    median_df_pred$lower),2),
                            "||",
                            "ELPD: ",
                            round(median_elpd_loo,2)))
p2
grid.arrange(p1,p2,nrow=1)

## calculation of residuals

## density function of SkewDoubleExponential distribution
## https://mc-stan.org/docs/functions-reference/skew-double-exponential-distribution.html
SkewDoubleExponential <- function(y, mu, sigma, tau) {
  result <- numeric(length(y))
  
  for (i in seq_along(y)) {
    if (y[i] < mu) {
      result[i] <- (2 * tau * (1 - tau) / sigma) * exp(-2 / sigma * (1 - tau) * (mu - y[i]))
    } else {
      result[i] <- (2 * tau * (1 - tau) / sigma) * exp(-2 / sigma * tau * (y[i] - mu))
    }
  }
  
  return(result)
}

sigma_summary <- median_stan_fit$summary("sigma")
df_ALD_residual <- data.frame(x = seq(from = -60, to = 40, length.out = 1000))
df_ALD_residual$density <- SkewDoubleExponential(y = df_ALD_residual$x,
                                                 mu = 0,
                                                 sigma = sigma_summary$mean,
                                                 tau = 0.5)
df_ALD_residual$loglik <- "ALD"

df_ALD_residual %>%
  ggplot(aes(x = x, y = density)) +
  geom_line()

## compile stan program for TPSC
TPSC_stan_file <- "../TPSC_t_MLR.stan"
TPSC_stan_mod <- cmdstan_model(TPSC_stan_file)

## MCMC modal regression
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

par.est <- TPSC_stan_fit$summary(c("alpha",paste0("beta[",1:P,"]")))
print(par.est)

## traceplot
pdf("traceplot_TPSC.pdf",height = 3,width = 8)
bayesplot::mcmc_trace(TPSC_stan_fit$draws(c("alpha",
                                            paste0("beta[",1:P,"]"),
                                            "w",
                                            "sigma",
                                            "delta")),
                      facet_args = list(ncol = 2))
dev.off()

## summary statistics
TPSC_stan_fit$summary(c("alpha",
                        paste0("beta[",1:P,"]"),
                        "w",
                        "sigma",
                        "delta")) %>%
  kbl(format = "latex", booktabs = TRUE, digits = 2) %>%
  kable_styling()

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

## modal fitting plot
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
  ggtitle("Likelihood: TPSC-Student-t (modal regression model)",
          subtitle = paste0("Coverage Rate:",
                            round(mean((y < TPSC_df_pred$upper) & 
                                         (y > TPSC_df_pred$lower)),2),
                            "||",
                            "Width: ",round(mean(TPSC_df_pred$upper -
                                                   TPSC_df_pred$lower),2),
                            "||",
                            "ELPD: ",
                            round(TPSC_elpd_loo,2)))
p3

## calculation of residuals
TPSC_student_t <- function(y,w,sigma,delta) {
  scale_LT <- sigma * sqrt(w/(1-w))
  scale_RT <- sigma * sqrt((1-w)/w)
  output <- numeric(length(y))
  for (i in seq_along(y)) {
    if (y[i] < 0) {
      output[i] <- w * 2 / scale_LT * dt(x = y[i]/scale_LT, df = delta)
    } else {
      output[i] <- (1 - w) * 2 / scale_RT * dt(x = y[i]/scale_RT, df = delta)
    }
  }
  return(output)
}

w_summary <- TPSC_stan_fit$summary("w")
delta_summary <- TPSC_stan_fit$summary("delta")
sigma_summary <- TPSC_stan_fit$summary("sigma")

df_TPSC_residual <- data.frame(x = seq(from = -60, to = 40, length.out = 1000))
df_TPSC_residual$density <- TPSC_student_t(y = df_TPSC_residual$x,
                                           w = w_summary$mean,
                                           sigma = sigma_summary$mean,
                                           delta = delta_summary$mean)
df_TPSC_residual$loglik <- "TPSC-Student-t"
df_TPSC_residual %>%
  ggplot(aes(x = x, y = density)) +
  geom_line()

## compile stan program for SNCP
SNCP_stan_file <- "../SN_MLR.stan"
SNCP_stan_mod <- cmdstan_model(SNCP_stan_file)

## MCMC robust mean regression
SNCP_stan_fit <- SNCP_stan_mod$sample(
  data = dat,
  seed = 100,
  chains = 4,
  parallel_chains = 4,
  refresh = 0,
  iter_warmup = 5000,
  save_warmup = FALSE,
  iter_sampling = 10000
)

par.est <- SNCP_stan_fit$summary(c("alpha",paste0("beta[",1:P,"]")))
print(par.est)

## traceplot
pdf("traceplot_SNCP.pdf",height = 3,width = 8)
bayesplot::mcmc_trace(SNCP_stan_fit$draws(c("alpha",
                                            paste0("beta[",1:P,"]"),
                                            "sigma",
                                            "gamma1")),
                      facet_args = list(ncol = 2))
dev.off()

SNCP_stan_fit$summary(c("alpha",
                        paste0("beta[",1:P,"]"),
                        "sigma",
                        "gamma1")) %>%
  kbl(format = "latex", booktabs = TRUE, digits = 2) %>%
  kable_styling()

## prediction interval
SNCP_stan_pred <- SNCP_stan_fit$draws(paste0("ystar[",1:N,"]"),
                                      format = "df")
SNCP_point_pred <- apply(SNCP_stan_pred, 2, median)[1:length(y)]
SNCP_hdi_df <- hdi(SNCP_stan_pred,credMass = 0.9)[,1:length(y)]
SNCP_df_pred <- data.frame(pred = SNCP_point_pred,
                           lower = SNCP_hdi_df[1,],
                           upper = SNCP_hdi_df[2,])

## model selection
SNCP_model_selection <- SNCP_stan_fit$loo()
SNCP_elpd_loo <- SNCP_model_selection$estimates[1,1]

## modal fitting plot
p4 <- df1 %>%
  ggplot() +
  geom_ribbon(data = SNCP_df_pred, 
              aes(x = df1$X, ymin = lower, ymax=upper),
              fill = "grey70") +
  geom_line(data = SNCP_df_pred, 
            aes(x = df1$X, y = pred),
            color = rgb(0.8,0,0,0.8)) +
  geom_point(aes(x = X, y = y)) +
  theme_bw() +
  ylim(-60,40) +
  xlab("X") +
  ylab("y") +
  ggtitle("Likelihood: SNCP (robust mean regression model)",
          subtitle = paste0("Coverage Rate:",
                            round(mean((y < SNCP_df_pred$upper) & 
                                         (y > SNCP_df_pred$lower)),2),
                            "||",
                            "Width: ",round(mean(SNCP_df_pred$upper -
                                                   SNCP_df_pred$lower),2),
                            "||",
                            "ELPD: ",
                            round(SNCP_elpd_loo,2)))
p4

## calculation of residuals

# The density function of skew normal with direct parameterization
# y: the value of one observation
# xi: the location parameter
# omega: the scale parameter
# alpha: the skewness parameter
dSNDP <- function(y,xi,omega,alpha) {
  p1 <- 2.0 / omega
  z <- (y - xi)/omega
  p2 <- dnorm(z)
  p3 <- pnorm(alpha*z)
  output <- p1 * p2 * p3
  return(output)
}

# The density function of skew normal with centred parameterization
# y: the value of one observation
# mu: the mean
# sigma: the standard deviation
# gamma1: the skewness parameter with range from -1 to 1
dSNCP <- function(y,mu,sigma,gamma1) {
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
  output <- dSNDP(y = y, xi = xi, omega = omega, alpha = alpha)
  return(output)
}

sigma_summary <- SNCP_stan_fit$summary("sigma")
gamma1_summary <- SNCP_stan_fit$summary("gamma1")
df_SNCP_residual <- data.frame(x = seq(from = -60, to = 40, length.out = 1000))
# the following calculation is temporary as the mode of SNCP is under calculation.
df_SNCP_residual$density <- dSNCP(y = df_SNCP_residual$x,
                                  mu = 0.0,
                                  sigma = sigma_summary$mean,
                                  gamma1 = gamma1_summary$mean)
# record the index of global mode of SNCP
index <- which.max(df_SNCP_residual$density)
# record the location shift value
loc_shift <- df_SNCP_residual$x[index]
# recalculate the density of SNCP with mode being 0
df_SNCP_residual$density <- dSNCP(y = df_SNCP_residual$x,
                                  mu = -loc_shift,
                                  sigma = sigma_summary$mean,
                                  gamma1 = gamma1_summary$mean)

df_SNCP_residual$loglik <- "SNCP"
df_SNCP_residual %>%
  ggplot(aes(x = x, y = density)) +
  geom_line()

## residual plot
df_residual_all <- rbind(df_error,
                         df_normal_residual,
                         df_SNCP_residual,
                         df_ALD_residual,
                         df_TPSC_residual)
unique(df_residual_all$loglik)
df_residual_all$loglik <- factor(df_residual_all$loglik,
                                 levels = c("True","Normal","SNCP","ALD","TPSC-Student-t"))

p_residual_all <- df_residual_all %>%
  mutate(Likelihood = loglik) %>%
  ggplot(aes(x = x , y = density)) +
  geom_line(aes(linetype = Likelihood, color = Likelihood)) +
  xlab("residual") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("left_skewed_residual.pdf", p_residual_all, width = 8, height = 4)


# plot of prediction ------------------------------------------------------

pall <- grid.arrange(p1,p4,p2,p3,ncol=1)
ggsave("left_skewed_data.pdf", pall, width = 8, height = 8)
