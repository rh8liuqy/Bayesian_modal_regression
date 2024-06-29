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

source("./MISC.R")


# simulate data -----------------------------------------------------------

y <- rTPSC(n = 100,w = 0.4,theta = 0,sigma = 1,delta = 5)


# complete data log-likelihood --------------------------------------------

dloglik <- function(theta,x,z,w = 0.4, sigma = 1, delta = 5){
  sigma1 <- sigma*sqrt(w/(1-w))
  sigma2 <- sigma*sqrt((1-w)/w)
  p1 <- w*2/sigma1*dtdist((x-theta)/sigma1,delta)*(x < theta)
  p1 <- p1^z
  p2 <- (1-w)*2/sigma2*dtdist((x-theta)/sigma2,delta)*(x >= theta)
  p2 <- p2^(1-z)
  output <- log(p1*p2)
  return(output)
}

# define the data augmented MCMC algorithm (reducible) --------------------

MCMC_algorithm <- function(theta_value, 
                           mcmc_size,
                           y) {
  z_output <- vector("list",mcmc_size)
  theta_output <- vector("list",mcmc_size)
  
  for (j in 1:mcmc_size) {

# update z ----------------------------------------------------------------
    
    z_value <- sapply(y, function(x){return(as.numeric(x < theta_value))})
    
# update theta ------------------------------------------------------------
    
    theta_value <- ESS(loglikelihood = function(theta) {
      output <- 0.0
      for (i in seq_along(y)) {
        output <- output + dloglik(theta = theta, 
                                   x = y[i], 
                                   z = z_value[i],
                                   w = 0.4, 
                                   sigma = 1, 
                                   delta = 5)
      }
      return(output)
    },
    f.current = theta_value,
    mu = 0,
    Sigma = 10^4)

# update all parameter values ---------------------------------------------

    z_output[[j]] <- z_value
    theta_output[[j]] <- theta_value
    
  }
  
  output <- list(z_output = z_output,
                 theta_output = theta_output)
  
  return(output)
  
}

MCMC_output <- MCMC_algorithm(theta_value = 8, 
                              mcmc_size = 2000,
                              y = y)

theta_output <- unlist(MCMC_output$theta_output)
df_theta <- data.frame(iteration = 1:2000,
                       theta = theta_output[1:2000])
p1 <- lattice::xyplot(theta ~ iteration, data = df_theta, 
                      type = "l",
                      main = "Initial value of theta is 8.")

MCMC_output <- MCMC_algorithm(theta_value = 1, 
                              mcmc_size = 2000,
                              y = y)

theta_output <- unlist(MCMC_output$theta_output)
df_theta <- data.frame(iteration = 1:2000,
                       theta = theta_output[1:2000])
p2 <- lattice::xyplot(theta ~ iteration, data = df_theta, 
                      type = "l",
                      main = "Initial value of theta is 1.")

MCMC_output <- MCMC_algorithm(theta_value = -5, 
                              mcmc_size = 2000,
                              y = y)

theta_output <- unlist(MCMC_output$theta_output)
df_theta <- data.frame(iteration = 1:2000,
                       theta = theta_output[1:2000])
p3 <- lattice::xyplot(theta ~ iteration, data = df_theta, 
                      type = "l",
                      main = "Initial value of theta is -5.")

pdf("./simulation_data_augmentation.pdf",width = 8.5, height = 8.5)
gridExtra::grid.arrange(p1,p2,p3,ncol = 1)
dev.off()