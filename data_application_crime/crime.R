## chapter 3.4, table 3.6
## https://plus.pearson.com/home
## https://users.stat.ufl.edu/~aa/social/data.html defines variable

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


# creation of scatter plot ------------------------------------------------

df <- read_xlsx("./crime.xlsx")

pdf("scatter_plot_US_crime_DC_include.pdf",width = 7,height = 7)
splom(df[,c(6,4,9,3)],main = NULL)
dev.off()

# Follow the suggestion from reviewer 1, remove District of Columbia in plot.

df <- df %>%
  filter(State != "District of Columbia")

pdf("scatter_plot_US_crime_DC_exclude.pdf",width = 7,height = 7)
splom(df[,c(6,4,9,3)],main = NULL)
dev.off()

data_application_crime <- function(option,remove_DC) {
  ## data import
  df1 <- read_xlsx("./crime.xlsx")
  if (remove_DC == TRUE) {
    df1 <- df1 %>%
      filter(State != "District of Columbia")
  }
  
  y <- df1$`murder rate`
  X <- as.matrix(df1[,c("college","poverty","metropolitan")])
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
  if (option == "normal" | option == "ALD") {
  df_post <- stan_fit$draws(variables = c("alpha",
                                          paste0("beta[",1:P,"]"),
                                          "sigma"), 
                            format = "df")
  } else if (option == "TPSC") {
    df_post <- stan_fit$draws(variables = c("alpha",
                                            paste0("beta[",1:P,"]"),
                                            "w",
                                            "sigma",
                                            "delta"), 
                              format = "df")
  }
  plot(mcmc_trace(df_post))
  
  ## LaTeX table in supplementary materials
  library(kableExtra)
  posterior_inference <- stan_fit$summary()
  print(posterior_inference %>%
          filter(! variable %in% c("lp__",
                                   paste0("ystar[",1:N,"]"),
                                   paste0("log_lik[",1:N,"]"))) %>%
          kbl(format = "latex", booktabs = TRUE, digits = 2) %>%
          kable_styling())
  
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
  
  return(list(par.est = par.est,
              elpd_loo = elpd_loo))
}

# Data analysis with DC ---------------------------------------------------

pdf("traceplot_normal_DC_included.pdf",height = 4,width = 8)
out1 <- data_application_crime("normal",remove_DC = FALSE)
dev.off()
pdf("traceplot_ALD_DC_included.pdf",height = 4,width = 8)
out2 <- data_application_crime("ALD",remove_DC = FALSE)
dev.off()
pdf("traceplot_TPSC_DC_included.pdf",height = 4,width = 8)
out3 <- data_application_crime("TPSC",remove_DC = FALSE)
dev.off()

out1
out2
out3

df_est1 <- bind_rows(out1$par.est[c("variable","mean","q5","q95")],
                     out2$par.est[c("variable","mean","q5","q95")],
                     out3$par.est[c("variable","mean","q5","q95")])
df_est1["Likelihood"] <- rep(c("normal","ALD","TPSC"),
                             each = 4)
df_est1["elpd_loo"] <- rep(c(out1$elpd_loo,
                             out2$elpd_loo,
                             out3$elpd_loo),
                           each = 4)
df_est1 <- df_est1[,c(5,6,1,2,3,4)]

# Data analysis without DC ------------------------------------------------

pdf("traceplot_normal_DC_excluded.pdf",height = 4,width = 8)
out4 <- data_application_crime("normal",remove_DC = TRUE)
dev.off()
pdf("traceplot_ALD_DC_excluded.pdf",height = 4,width = 8)
out5 <- data_application_crime("ALD",remove_DC = TRUE)
dev.off()
pdf("traceplot_TPSC_DC_excluded.pdf",height = 4,width = 8)
out6 <- data_application_crime("TPSC",remove_DC = TRUE)
dev.off()

out4
out5
out6

df_est2 <- bind_rows(out4$par.est[c("variable","mean","q5","q95")],
                     out5$par.est[c("variable","mean","q5","q95")],
                     out6$par.est[c("variable","mean","q5","q95")])
df_est2["Likelihood"] <- rep(c("normal","ALD","TPSC"),
                             each = 4)
df_est2["elpd_loo"] <- rep(c(out4$elpd_loo,
                             out5$elpd_loo,
                             out6$elpd_loo),
                           each = 4)
df_est2 <- df_est2[,c(5,6,1,2,3,4)]
df_est <- rbind(df_est1,df_est2)
df_est$DC <- rep(c("Include D.C.",
                   "Exclude D.C."),
                 each = 12)
df_est <- df_est[,c(7,1:6)]

kbl(df_est,
    digits = 4,
    format = "latex",
    booktabs = TRUE,
    caption = "this is caption",
    label = "this is label") %>%
  collapse_rows(columns = 1:3, valign = "middle")

kbl(df_est,
    digits = 4,
    booktabs = TRUE,
    caption = "this is caption",
    label = "this is label") %>%
  collapse_rows(columns = 1:3, valign = "middle") %>%
  kable_styling()
