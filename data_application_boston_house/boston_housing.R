rm(list = ls())
library(cmdstanr)
library(posterior)
library(bayesplot)
library(MASS)
library(tidyverse)
library(HDInterval)
library(kableExtra)

set.seed(100)

data_application <- function(option) {
  df1 <- MASS::Boston
  
  ## data list for stan
  N <- nrow(df1)
  y <- df1$medv
  X <- df1 %>%
    dplyr::select(-medv)
  X <- as.matrix(X)
  P <- ncol(X)
  
  dat <- list(N = N,
              X = X,
              y = y,
              P = P)
  
  if (option == "normal") {
    ## compile stan program for normal
    stan_file <- "../normal_MLR.stan"
    title.txt <- "Likelihood: Normal (mean regression model)"
  }
  else if(option == "ALD") {
    stan_file <- "../quantile_MLR.stan"
    title.txt <- "Likelihood: ALD (median regression model)"
  }
  else if(option == "TPSC") {
    stan_file <- "../TPSC_t_MLR.stan"
    title.txt <- "Likelihood: TPSC-Student-t (modal regression model)"
  }
  else if(option == "FG") {
    stan_file <- "../FG_MLR.stan"
    title.txt <- "Likelihood: FG (modal regression model)"
  }
  else if(option == "SN") {
    stan_file <- "../SN_MLR.stan"
    title.txt <- "Likelihood: SN_CP (mean regression model)"
  }
  
  stan_mod <- cmdstan_model(stan_file)
  
  ## load initial value
  df_initial <- read.csv("./initial_value.csv")
  list_initial <- vector(mode = "list", length = 2)
  names(list_initial) <- c("alpha","beta")
  list_initial[[1]] <- df_initial$median[1]
  list_initial[[2]] <- df_initial$median[2:nrow(df_initial)]
  if (option == "SN") {
    list_initial[[3]] <- 1.0
    list_initial[[4]] <- 0.0
    names(list_initial)[3:4] <- c("sigma","gamma1")
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
    show_exceptions = FALSE,
    init = function() {return(list_initial)}
  )
  
  ## parameter estimation 
  par.est <- stan_fit$summary(c("alpha",paste0("beta[",1:P,"]")))
  write.csv(x = par.est, 
            file = paste0(option,"_inference.csv"),
            row.names = FALSE)
  
  ## traceplots and summary statistics
  if (option == "normal" | option == "ALD") {
    df_post <- stan_fit$draws(variables = c("alpha",
                                            paste0("beta[",1:P,"]"),
                                            "sigma"), 
                              format = "df") 
    posterior_inference <- stan_fit$summary(c("alpha",
                                              paste0("beta[",1:P,"]"),
                                              "sigma"))
  } else if (option == "SN"){
    df_post <- stan_fit$draws(variables = c("alpha",
                                            paste0("beta[",1:P,"]"),
                                            "sigma",
                                            "gamma1"), 
                              format = "df")
    posterior_inference <- stan_fit$summary(c("alpha",
                                              paste0("beta[",1:P,"]"),
                                              "sigma",
                                              "gamma1"))
  } else if (option == "TPSC"){
    df_post <- stan_fit$draws(variables = c("alpha",
                                            paste0("beta[",1:P,"]"),
                                            "w",
                                            "sigma",
                                            "delta"), 
                              format = "df") 
    posterior_inference <- stan_fit$summary(c("alpha",
                                              paste0("beta[",1:P,"]"),
                                              "w",
                                              "sigma",
                                              "delta"))
  }
  
  ## print MCMC trace
  plot(mcmc_trace(df_post, facet_args = list(ncol = 2)))
  
  ## print summary statistics
  print(posterior_inference %>%
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
  
  ## print likelihood type
  print(title.txt)
  
  output <- list(coverage = round(mean((y < df_pred$upper) & (y > df_pred$lower)), 4),
                 width = round(mean(df_pred$upper - df_pred$lower),2),
                 ELPD = round(elpd_loo,2))
  ## Coverage Rate
  print(paste0("Coverage Rate: ", output$coverage))
  ## Width
  print(paste0("Width: ", output$width))
  ## ELPD
  print(paste0("ELPD: ", output$ELPD))
  return(output)
}

pdf("traceplot_normal.pdf",height = 11,width = 8.5)
output_normal <- data_application("normal")
dev.off()
pdf("traceplot_SN.pdf",height = 11,width = 8.5)
output_SN <- data_application("SN")
dev.off()
pdf("traceplot_ALD.pdf",height = 11,width = 8.5)
output_ALD <- data_application("ALD")
dev.off()
pdf("traceplot_TPSC.pdf",height = 11,width = 8.5)
output_TPSC <- data_application("TPSC")
dev.off()


# Table for Coverage Rate/Width/ELPD --------------------------------------
df_tab <- as.data.frame(matrix(0, nrow = 4, ncol = 3))
rownames(df_tab) <- c("Normal","SNCP","ALD","TPSC-Student-t")
colnames(df_tab) <- c("Coverage Rate (%)", "Width", "ELPD")
df_tab[1,] <- unlist(output_normal)
df_tab[2,] <- unlist(output_SN)
df_tab[3,] <- unlist(output_ALD)
df_tab[4,] <- unlist(output_TPSC)
df_tab[,1] <- 100 * df_tab[,1]
df_tab %>%
  kbl(booktabs = TRUE, format = "latex") %>%
  kable_classic()

# Inference Table for betas -----------------------------------------------

## inference table - normal
output <- read.csv("./normal_inference.csv")
output_normal <- output %>%
  mutate(across(2:10, \(x) round(x,2))) %>%
  select(variable,mean,q5,q95)

## inference table - ALD
output <- read.csv("./ALD_inference.csv")
output_ALD <- output %>%
  mutate(across(2:10, \(x) round(x,2))) %>%
  select(mean,q5,q95)

## inference table - SN
output <- read.csv("./SN_inference.csv")
output_SN <- output %>%
  mutate(across(2:10, \(x) round(x,2))) %>%
  select(mean,q5,q95)

## inference table - TPSC
output <- read.csv("./TPSC_inference.csv")
output_TPSC <- output %>%
  mutate(across(2:10, \(x) round(x,2))) %>%
  select(mean,q5,q95)

df_tab <- cbind(output_normal, output_SN, output_ALD, output_TPSC)
colnames(df_tab) <- NULL

df_tab %>%
  kbl(format = "latex", booktabs = TRUE) %>%
  kable_classic() %>%
  add_header_above(c("Parameter" = 1, 
                     "Mean" = 1, "q5" = 1, "q95" = 1,
                     "Mean" = 1, "q5" = 1, "q95" = 1,
                     "Mean" = 1, "q5" = 1, "q95" = 1,
                     "Mean" = 1, "q5" = 1, "q95" = 1)) %>%
  add_header_above(c(" ", "Normal" = 3, "SNCP" = 3, "ALD" = 3, "TPSC-Student-t" = 3))
