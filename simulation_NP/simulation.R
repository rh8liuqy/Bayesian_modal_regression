rm(list = ls())

library(kableExtra)

source("./Model_Weixin_Yao.R")

set.seed(100)

# approximate the mode of SN(0,1,5) ---------------------------------------
# SN(0,1,5) stands for the SN distribution (direct parameterization) with location 0, scale 1 and skewness 5.

x <- seq(-5,5,length.out = 1000)
value <- sn::dsn(x = x,xi = 0, omega = 1,alpha = 5)

SN_mode <- x[which.max(value)]

# determining the optimal bandwidth ---------------------------------------

simu_MODLR_h <- function(seed,h) {
  set.seed(seed)
  n <- 30
  beta_true <- rep(1,3)
  
  noise <- sn::rsn(n = n, xi = 0, omega = 1, alpha = 5) - SN_mode
  
  X <- matrix(data = rnorm(n*2), ncol = 2)
  X <- cbind(1,X)
  eta <- as.numeric(X %*% beta_true)
  y <- eta + noise
  
  output <- MODLR(y = y, X = X, beta_init = rnorm(ncol(X)), h = h, iter_max = 1000, tol = 1e-8)

  output <- matrix(output,nrow = 1)
  return(output)
}

# define the vector and bandwidth

h_vec <- seq(0.3,1,length.out = 100)

# find the optimal bandwidth across 300 Monte Carlo replicates

bias_output <- matrix(data = NA, nrow = 100, ncol = 3)

i <- 1
for (hi in h_vec) {
  simu_MODLR_output <- parallel::mclapply(X = 1:300,FUN = simu_MODLR_h, h = hi, mc.cores = 8)
  simu_MODLR_output <- do.call(rbind,simu_MODLR_output)
  bias_output[i,] <- colMeans(simu_MODLR_output) - rep(1,3)
  i <- i + 1
}

index_h <- which.min(apply(bias_output, MARGIN = 1, function(x) {sum(x^2)}))

# h_opt is the optimal bandwidth

h_opt <- h_vec[index_h]

# repeat the simulation with the optimal bandwidth ------------------------

simu_MODLR <- function(seed,h) {
  set.seed(seed)
  n <- 30
  beta_true <- rep(1,3)
  
  noise <- sn::rsn(n = n, xi = 0, omega = 1, alpha = 5) - SN_mode
  
  X <- matrix(data = rnorm(n*2), ncol = 2)
  X <- cbind(1,X)
  eta <- as.numeric(X %*% beta_true)
  y <- eta + noise
  
  # point estimation
  point_est <- MODLR(y = y, X = X, beta_init = beta_true, h = h, iter_max = 1000, tol = 1e-8)
  
  # define residuals
  MODLR_residuals <- as.numeric(y - X %*% point_est)
  
  # bootstrap for 500 times
  B <- 500
  output_B <- matrix(data = NA, nrow = B, ncol = ncol(X))
  i <- 1
  while (i <= B) {
    bootstrap_index <- sample(1:n,size = n, replace = TRUE)
    
    yB <- y[bootstrap_index]
    XB <- X[bootstrap_index,]
    
    # point estimation may fail because of numerical issues
    
    point_est_B <- tryCatch({
      MODLR(y = yB, X = XB, beta_init = point_est, h = h, iter_max = 1000, tol = 1e-8)
    }, error = function(e) {return(NULL)})
    
    # if point estimation is successful then move on
    
    if (! is.null(point_est_B)) {
      output_B[i,] <- point_est_B
      i <- i + 1
    }
  }
  
  # lower bound of betas
  lb <- apply(output_B, MARGIN = 2, function(x) {quantile(x,probs = 0.05)})
  # upper bound of betas
  ub <- apply(output_B, MARGIN = 2, function(x) {quantile(x,probs = 0.95)})
  
  output <- data.frame(beta0_PE = point_est[1],
                       beta1_PE = point_est[2],
                       beta2_PE = point_est[3],
                       beta0_CO = lb[1] < 1 & 1 < ub[1],
                       beta1_CO = lb[2] < 1 & 1 < ub[2],
                       beta2_CO = lb[3] < 1 & 1 < ub[3],
                       beta0_CL = ub[1] - lb[1],
                       beta1_CL = ub[2] - lb[2],
                       beta2_CL = ub[3] - lb[3])
  return(output)
}

simu_MODLR_output <- parallel::mclapply(X = 1:300,FUN = simu_MODLR, h = h_opt, mc.cores = 8)
simu_MODLR_output <- do.call(rbind,simu_MODLR_output)

# Repeat the same experiment with GUD modal regression --------------------

library(cmdstanr)

## compile stan program for TPSC
TPSC_stan_file <- "../TPSC_t_MLR.stan"
TPSC_stan_mod <- cmdstan_model(TPSC_stan_file)

simu_GUD <- function(seed) {
  set.seed(seed)
  n <- 30
  beta_true <- rep(1,3)
  
  noise <- sn::rsn(n = n, xi = 0, omega = 1, alpha = 5) - SN_mode
  
  X <- matrix(data = rnorm(n*2), ncol = 2)
  X <- cbind(1,X)
  eta <- as.numeric(X %*% beta_true)
  y <- eta + noise
  
  # remove the intercept
  X <- X[,-1]
  
  dat <- list(N = n,
              X = X,
              P = 2,
              y = y)
  
  TPSC_stan_fit <- TPSC_stan_mod$sample(
    data = dat,
    seed = 100,
    chains = 1,
    parallel_chains = 1,
    refresh = 0,
    iter_warmup = 5000,
    save_warmup = FALSE,
    iter_sampling = 10000
  )
  
  par.est <- TPSC_stan_fit$summary(c("alpha","beta"))
  
  df_alpha <- TPSC_stan_fit$draws("alpha",format = "df")
  df_alpha <- as.data.frame(df_alpha)[,1]
  alpha_ci <- HDInterval::hdi(df_alpha,credMass = 0.90)
  
  df_beta1 <- TPSC_stan_fit$draws("beta[1]",format = "df")
  df_beta1 <- as.data.frame(df_beta1)[,1]
  beta1_ci <- HDInterval::hdi(df_beta1,credMass = 0.90)
  
  df_beta2 <- TPSC_stan_fit$draws("beta[2]",format = "df")
  df_beta2 <- as.data.frame(df_beta2)[,1]
  beta2_ci <- HDInterval::hdi(df_beta2,credMass = 0.90)
  
  output <- data.frame(beta0_PE = par.est$mean[1],
                       beta1_PE = par.est$mean[2],
                       beta2_PE = par.est$mean[3],
                       beta0_CO = alpha_ci[1] < 1 & 1 < alpha_ci[2],
                       beta1_CO = beta2_ci[1] < 1 & 1 < beta1_ci[2],
                       beta2_CO = beta2_ci[1] < 1 & 1 < beta2_ci[2],
                       beta0_CL = alpha_ci[2] - alpha_ci[1],
                       beta1_CL = beta1_ci[2] - beta1_ci[1],
                       beta2_CL = beta2_ci[2] - beta2_ci[1])
  return(output)
}

simu_GUD_output <- parallel::mclapply(X = 1:300,FUN = simu_GUD, mc.cores = 8)
simu_GUD_output <- do.call(rbind,simu_GUD_output)

# point estimation and coverage rate - nonparametric modal regression

avg_MODLR <- colMeans(simu_MODLR_output)

# Monte-Carlo standard deviation - nonparametric modal regression

sd_MODLR <- apply(simu_MODLR_output[,c(1:3,7:9)], MARGIN = 2, function(x){sd(x)})

# table for MODLR 

df_MODLR <- data.frame(Regression.Model = rep("MODLR",3),
                       Parameter = c("beta0","beta1","beta2"),
                       Average.Point.Estimation = c(paste0(round(avg_MODLR[1],2),"(",round(sd_MODLR[1],2),")"),
                                                    paste0(round(avg_MODLR[2],2),"(",round(sd_MODLR[2],2),")"),
                                                    paste0(round(avg_MODLR[3],2),"(",round(sd_MODLR[3],2),")")),
                       Coverage.Rate = round(avg_MODLR[4:6]*100,2),
                       Width.CI = c(paste0(round(avg_MODLR[7],2),"(",round(sd_MODLR[4],2),")"),
                                    paste0(round(avg_MODLR[8],2),"(",round(sd_MODLR[5],2),")"),
                                    paste0(round(avg_MODLR[9],2),"(",round(sd_MODLR[6],2),")")))

rownames(df_MODLR) <- NULL

# point estimation and coverage rate - GUD regression

avg_GUD <- colMeans(simu_GUD_output)

# Monte-Carlo standard deviation - GUD regression

sd_GUD <- apply(simu_GUD_output[,c(1:3,7:9)], MARGIN = 2, function(x){sd(x)})

# table for GUD

df_GUD <- data.frame(Regression.Model = rep("Bayesian Modal Regression",3),
                     Parameter = c("beta0","beta1","beta2"),
                     Average.Point.Estimation = c(paste0(round(avg_GUD[1],2),"(",round(sd_GUD[1],2),")"),
                                                  paste0(round(avg_GUD[2],2),"(",round(sd_GUD[2],2),")"),
                                                  paste0(round(avg_GUD[3],2),"(",round(sd_GUD[3],2),")")),
                     Coverage.Rate = round(avg_GUD[4:6]*100,2),
                     Width.CI = c(paste0(round(avg_GUD[7],2),"(",round(sd_GUD[4],2),")"),
                                  paste0(round(avg_GUD[8],2),"(",round(sd_GUD[5],2),")"),
                                  paste0(round(avg_GUD[9],2),"(",round(sd_GUD[6],2),")")))

rownames(df_GUD) <- NULL

# combine table for output
rbind(df_GUD,df_MODLR) %>%
  kbl(format = "latex",booktabs = TRUE) %>%
  collapse_rows(columns = c(1,2)) %>%
  kable_styling()
