## Data from Yu, K. and Moyeed, R. A. (2001), "Bayesian Quantile Regression," Statistics and Probability Letters, 54, 437–447.
## Model evaluation: Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC.
## elpd_loo: leave-one-out expected log predictive density
## p_loo: effective number of parameters
## looic: -2*elpd_loo

rm(list = ls())
library(cmdstanr)
library(posterior)
library(bayesplot)
library(coda)
library(MASS)
library(tidyverse)
library(readxl)
library(gridExtra)
library(latex2exp)
library(parallel)
library(HDInterval)
color_scheme_set("brightblue")

data_application_serum <- function(option) {
  ## data import
  df1 <- readxl::read_xlsx("./IgG.xlsx")
  df1$Age2 <- (df1$Age)^2
  
  y <- df1$IgG
  X <- as.matrix(df1[,c("Age","Age2")])
  N <- nrow(X)
  P <- ncol(X)
  
  ## data list for stan
  
  dat <- list(N = N,
              P = P,
              X = X,
              y = y)
  if (option == "normal") {
    stan_file <- "../normal_MLR.stan"
    stan_mod <- cmdstan_model(stan_file)
    title.txt <- "Likelihood: Normal (mean regression model)"
  }
  else if(option == "ALD") {
    stan_file <- "../quantile_MLR.stan"
    stan_mod <- cmdstan_model(stan_file)
    title.txt <- "Likelihood: ALD (median regression model)"
  }
  else if(option == "FG") {
    stan_file <- "../FG_MLR.stan"
    stan_mod <- cmdstan_model(stan_file)
    title.txt <- "Likelihood: FG (modal regression model)"
  }
  else if(option == "logNM") {
    stan_file <- "../logNM_MLR.stan"
    stan_mod <- cmdstan_model(stan_file)
    title.txt <- "Likelihood: logNM (modal regression model)"
  }
  else if(option == "TPSC") {
    stan_file <- "../TPSC_t_MLR.stan"
    stan_mod <- cmdstan_model(stan_file)
    title.txt <- "Likelihood: TPSC student t (modal regression model)"
  }
  else if(option == "DTP") {
    stan_file <- "../DTP_MLR.stan"
    stan_mod <- cmdstan_model(stan_file)
    title.txt <- "Likelihood: DTP student t (modal regression model)"
  }
  
  ## MCMC
  stan_fit <- stan_mod$sample(
    data = dat,
    seed = 100,
    chains = 4,
    parallel_chains = 4,
    refresh = 0,
    iter_warmup = 5000,
    save_warmup = FALSE,
    iter_sampling = 10000,
  )
  
  ## parameter estimation 
  par.est <- stan_fit$summary(c("alpha",paste0("beta[",1:P,"]")))
  print(par.est)
  df_post <- stan_fit$draws(variables = c("alpha",paste0("beta[",1:P,"]")), 
                            format = "df")
  plot(mcmc_trace(df_post))
  
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
  
  p1 <- df1 %>%
    ggplot() +
    geom_ribbon(data = df_pred, 
                aes(x = df1$Age, ymin = lower, ymax=upper),
                fill = "grey70") +
    geom_line(data = df_pred, 
              aes(x = df1$Age, y = pred),
              color = rgb(0.8,0,0,0.8)) +
    geom_point(aes(x = Age, y = IgG)) +
    theme_bw() +
    ylim(-2.5,15) +
    xlab("Age") +
    ylab("IgG") +
    ggtitle(title.txt,
            subtitle = paste0("Coverage Rate:",
                              round(mean((y < df_pred$upper) & 
                                           (y > df_pred$lower)),4),
                              "||",
                              "height: ",round(mean(df_pred$upper -
                                                      df_pred$lower),4),
                              "||",
                              "ELPD: ",
                              round(elpd_loo,4)))
  return(list(par.est = par.est,
              elpd_loo = elpd_loo,
              figure = p1))
}

pdf("traceplot_normal.pdf",height = 4,width = 8)
out1 <- data_application_serum("normal")
dev.off()
pdf("traceplot_ALD.pdf",height = 4,width = 8)
out2 <- data_application_serum("ALD")
dev.off()
pdf("traceplot_TPSC.pdf",height = 4,width = 8)
out3 <- data_application_serum("FG")
dev.off()
#out4 <- data_application_serum("logNM")
#out5 <- data_application_serum("TPSC")
#out6 <- data_application_serum("DTP")

out1
out2
out3

pall <- grid.arrange(out1$figure,
                     out2$figure,
                     out3$figure,
                     ncol = 2)
ggsave("serum.pdf",pall,width = 10,height = 10)

df_est <- bind_rows(out1$par.est[c("variable","mean","q5","q95")],
                    out2$par.est[c("variable","mean","q5","q95")],
                    out3$par.est[c("variable","mean","q5","q95")]
                    )
df_est["Likelihood"] <- rep(c("normal","ALD",
                              "FG"),
                            each = 3)
df_est$elpd_loo <- rep(c(out1$elpd_loo,
                         out2$elpd_loo,
                         out3$elpd_loo),
                       each = 3)
df_est <- df_est[,c(5,6,1,2,3,4)]

library(kableExtra)
kbl(df_est,
    digits = 4,
    format = "latex",
    booktabs = TRUE,
    caption = "this is caption",
    label = "this is label") %>%
  collapse_rows(columns = 1:2, valign = "middle")
