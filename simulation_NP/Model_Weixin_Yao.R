# Modal EM algorithm for nonparametric model regression by Weixin  --------
# A New Regression Model: Modal Linear Regression
# doi: 10.1111/sjos.12054

# The function for calculate the normal kernel with given bandwidth.
# h: bandwidth.

phi_h <- function(x,h) {
  output <- dnorm(x/h)/h
  return(output)
}

# The function for normalize a vector of proportional probabilities whose sum is not 1.
# probs: The vector of proportional probability.

norm_prob <- function(probs) {
  sum_probs <- sum(probs)
  output <- probs / sum_probs
  return(output)
}

# y: response variable y.
# X: the design matrix.
# beta_init: the initial value of betas.
# h: the bandwidth.
# iter_max: the maximum number of iterations allowed.
# tol: the tolerance level for convergence.
MODLR <- function(y,X,beta_init,h,iter_max,tol) {
  if (ncol(X) != length(beta_init)) {
    stop("Check the dimension of design matrix and initial values of beta.")
  }
  
  # define the sample size
  n <- length(y)
  
  # initialize the beta_value
  beta_value <- beta_init
  
  # the transpose of the design matrix
  Xt <- t(X)
  
  # initial the iteration number
  iter_num <- 1
  
  # initial the difference of beta between each iteration
  
  tol_M <- 1e3
  
  while (iter_num <= iter_max & tol_M > tol) {

# E step ------------------------------------------------------------------

    eta_value <- as.numeric(X %*% beta_value)
    pi_vec <- phi_h(y - eta_value,h)
    pi_vec <- norm_prob(pi_vec)

# M step ------------------------------------------------------------------

    W <- diag(x = pi_vec)
    XtW <- Xt %*% W
    XtWX <- XtW %*% X
    XtWy <- XtW %*% y
    beta_value_old <- beta_value
    beta_value <- solve(a = XtWX, b = XtWy)
    
    # define the change in beta values
    
    tol_M <- norm(beta_value_old - beta_value, type = "2")
    
    print(paste0("iter_num: ", iter_num, " || tol_M: ", tol_M))
    
    # recount the number of iterations
    
    iter_num <- iter_num + 1
    
  }
  output <- beta_value
  return(output)
}


# Bandwidth selection -----------------------------------------------------

# the kernel function for nu = 0

K0 <- function(x) {
  return(dnorm(x))
}

# the kernel function for nu = 3

K3 <- function(x) {
  return(dnorm(x)*(-x)*(x^2-3))
}

# The function for finding the best bandwidth.
# y: response variable y.
# X: the design matrix.
# residuals_value: The value of residuals.
hopt_func <- function(y,X,residuals_value) {

# define the sample size, n.
  
  n <- length(y)
  
# define the number of columns in the design matrix X
  
  p <- ncol(X)

# estimation the mode of least square residuals ---------------------------
  
  residuals_KDE <- density(x = residuals_value,
                           bw = "SJ",
                           kernel = "gaussian")
  
  # the bandwidth h
  
  h <- residuals_KDE$bw
  
  index_mode <- which.max(residuals_KDE$y)
  
  # m_hat is the estimated mode of residuals
  
  m_hat <- residuals_KDE$x[index_mode]

# estimation of K_hat and L_hat -------------------------------------------
  
  K_hat_mat <- matrix(NA,nrow = n,ncol = p)
  L_hat_mat <-matrix(0,nrow = p,ncol = p)
  
  for (i in 1:n) {
    
    K_hat_mat[i,] <- K3(x = (residuals_value[i] - m_hat)/h ) / (h^(3+1)) * as.numeric(X[i,])
    
    L_hat_mat <- L_hat_mat + K0(x = (residuals_value[i] - m_hat)/h ) / (h^(0+1)) * (matrix(X[i,],nrow = p) %*% matrix(X[i,],ncol = p))
    
  }
  
  K_hat <- colMeans(K_hat_mat)
  L_hat <- L_hat_mat/n

# define of the constant nu2 ----------------------------------------------
  
  # nu2 = int_{-Inf}^{Inf} t^2 (dnorm(t))^2 dt
  nu2 <- 1/(4*sqrt(pi))


# the calculation of the optimal bandwidth --------------------------------
  
  hopt <- ((3 * nu2 * (p+1))/(t(K_hat) %*% solve(L_hat, K_hat)))^(1/7) * n^(-1/7)
  
  output <- as.numeric(hopt)
  
  return(output)

}

# simulation to verify correctness ----------------------------------------

# set.seed(300)
# 
# simu <- function(temp) {
#   n <- 200
#   noise <- rt(n = n, df = 4)
#   X <- matrix(data = rnorm(n*4), ncol = 4)
#   beta_true <- rep(1,4)
#   eta_true <- as.numeric(X %*% beta_true)
#   y <- eta_true + noise
#   
#   beta_init <- rnorm(4)
#   
#   h1 <- hopt_func(y,X,residuals(lm(y ~ X - 1)))
# 
#   est1 <- MODLR(y = y, X = X, beta_init = beta_init, h = h1, iter_max = 1000, tol = 1e-8)
# 
#   residuals1 <- as.numeric(y - X %*% est1)
# 
#   h2 <- hopt_func(y,X,residuals1)
# 
#   est2 <- MODLR(y = y, X = X, beta_init = est1, h = h2, iter_max = 1000, tol = 1e-8)
#   
#   output <- matrix(est2,nrow = 1)
#   return(output)
# }
# 
# simu_output <- parallel::mclapply(1:1000,simu,mc.cores = 8)
# simu_output <- do.call(rbind,simu_output)
# colMeans(simu_output)
