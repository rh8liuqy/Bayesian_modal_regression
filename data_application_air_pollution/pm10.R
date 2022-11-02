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

set.seed(100)

data_application_NO2 <- function(option){
  ## data import
  df1 <- read.csv("./pm10.csv")
  
  ## data list for stan
  N <- nrow(df1)
  y <- df1$PM10
  X <- as.matrix(df1$wind_speed)
  P <- ncol(X)
  
  dat <- list(N = N,
              X = X,
              y = y,
              P = P)
  if (option == "normal") {
    ## compile stan program for normal
    stan_file <- "../normal_MLR.stan"
    stan_mod <- cmdstan_model(stan_file)
    title.txt <- "Likelihood: Normal (mean regression model)"
  }
  else if(option == "ALD") {
    stan_file <- "../quantile_MLR.stan"
    stan_mod <- cmdstan_model(stan_file)
    title.txt <- "Likelihood: ALD (median regression model)"
  }
  else if(option == "TPSC") {
    stan_file <- "../TPSC_t_MLR.stan"
    stan_mod <- cmdstan_model(stan_file)
    title.txt <- "Likelihood: TPSC student t (modal regression model)"
  }
  else if(option == "FG") {
    stan_file <- "../FG_MLR.stan"
    stan_mod <- cmdstan_model(stan_file)
    title.txt <- "Likelihood: FG (modal regression model)"
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
    # init = function(chain_id) {
    #     list( w1=w1,
    #           scale1=scale1,
    #           scale2=scale2,
    #           rho=1.3,
    #           alpha=1.4,
    #           b0=0.0,
    #           eta=rep(0,N))},
  )
  
  ## parameter estimation 
  par.est <- stan_fit$summary(c("alpha",paste0("beta[",1:P,"]")))
  print(par.est)
  
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
                aes(x = df1$wind_speed, ymin = lower, ymax=upper),
                fill = "grey",
                alpha = 0.4) +
    geom_line(data = df_pred, 
              aes(x = df1$wind_speed, y = pred),
              color = rgb(0.8,0,0,0.8)) +
    geom_point(aes(x = wind_speed, y = PM10)) +
    geom_hline(yintercept = 0,
               color = rgb(0,0,0.8,0.8),
               linetype = 2) +
    theme_bw() +
    xlab("wind speed") +
    ylab("PM10") +
    ylim(-50,300) +
    ggtitle(title.txt,
            subtitle = paste0("Coverage Rate:",
                              round(mean((y < df_pred$upper) & 
                                           (y > df_pred$lower)),4),
                              "||",
                              "height: ",round(mean(df_pred$upper -
                                                      df_pred$lower),4),
                              "||",
                              "elpd_loo: ",
                              round(elpd_loo,4)))
  return(return(list(par.est = par.est,
                     elpd_loo = elpd_loo,
                     figure = p1)))
}

out1 <- data_application_NO2("normal")
out2 <- data_application_NO2("ALD")
out3 <- data_application_NO2("TPSC")
out4 <- data_application_NO2("FG")

pall <- grid.arrange(out1$figure,out2$figure,out3$figure,nrow = 3)
ggsave("pm10.pdf",pall,width = 8,height = 8)
