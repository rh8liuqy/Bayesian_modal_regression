#' The Density Function of the Student-t Distribution
#'
#' @param x vector of quantiles.
#' @param delta vector of degree of freedom parameters.
#'
#' @return vector of densities.
#' @noRd
dtdist <- function(x,delta) {
  p1 <- gamma((delta+1)/2)/gamma(delta/2)
  p2 <- 1/sqrt(delta*pi)
  p3 <- (1+1/delta*(x)^2)^(-(delta+1)/2)
  return(p1*p2*p3)
}

#' The Random Number Generation of the Student-t Distribution
#'
#' @param n vector of quantiles.
#' @param delta vector of degree of freedom parameters.
#'
#' @return vector of random deviates.
#' @noRd
rtdist <- function(n,delta) {
  U <- rgamma(n = n, shape = delta/2, rate = delta/2)
  X <- rnorm(n = n)
  output <- U^(-1/2) * X
  return(output)
}

#' The TPSC-Student-t Distribution
#' @param x vector of quantiles.
#' @param n number of observations.
#' @param w vector of weight parameters.
#' @param theta vector of the location parameters.
#' @param sigma vector of the scale parameters.
#' @param delta the degree of freedom.
#' @name dTPSC
#' @details
#' The TPSC-Student-t distribution has density
#' \deqn{f_{\mathrm{TPSC}}(y \mid w, \theta, \sigma, \delta)=w f_{\mathrm{LT}}\left(y \mid \theta, \sigma \sqrt{\frac{w}{1-w}}, \delta\right)+(1-w) f_{\mathrm{RT}}\left(y \mid \theta, \sigma \sqrt{\frac{1-w}{w}}, \delta\right),}
#' where \deqn{f_{\mathrm{LT}}(y \mid \theta, \sigma, \delta)=\frac{2}{\sigma} f\left(\left.\frac{y-\theta}{\sigma} \right\rvert\, \delta\right) \mathbb{I}(y<\theta),} and \deqn{f_{\mathrm{RT}}(y \mid \theta, \sigma, \delta)=\frac{2}{\sigma} f\left(\left.\frac{y-\theta}{\sigma} \right\rvert\, \delta\right) \mathbb{I}(y \geq \theta).}
#'
#' @return `dTPSC` gives the density. `rTPSC` generates random deviates.
#' @example /inst/examples/TPSC_example.R
NULL
#> NULL

#' @rdname dTPSC
#' @export
dTPSC <- function(x,w,theta,sigma,delta){
  sigma1 <- sigma*sqrt(w/(1-w))
  sigma2 <- sigma*sqrt((1-w)/w)
  p1 <- w*2/sigma1*dtdist((x-theta)/sigma1,delta)*(x<theta)
  p2 <- (1-w)*2/sigma2*dtdist((x-theta)/sigma2,delta)*(x>=theta)
  return(p1+p2)
}

#' @rdname dTPSC
#' @export
rTPSC <- function(n,w,theta,sigma,delta){
  X1 <- -abs(rtdist(n = n,delta = delta)) * sigma * sqrt(w/(1-w))
  X2 <- abs(rtdist(n = n,delta = delta)) * sigma * sqrt((1-w)/w)
  Z <- rbinom(n = n, size = 1, prob = w)
  output <- X1*Z + X2*(1-Z) + theta
  return(output)
}

# The Elliptical Slice Sampling -------------------------------------------

ESS <- function(loglikelihood,
                f.current,
                mu,
                Sigma) {
  p <- length(f.current)
  nu <- mvrnorm(n = 1,
                mu = mu,
                Sigma = Sigma)
  nu <- as.numeric(nu)
  u <- runif(n = 1,
             min = 0,
             max = 1)
  logy <- loglikelihood(f.current) + log(u)
  theta <- runif(n = 1,
                 min = 0,
                 max = 2*pi)
  theta.min <- theta - 2*pi
  theta.max <- theta
  condition <- TRUE
  while (condition) {
    f.new <- (f.current-mu)*cos(theta) + (nu-mu)*sin(theta) + mu
    if (loglikelihood(f.new) > logy) {
      condition <- FALSE
    }
    else {
      if (theta < 0) {
        theta.min <- theta
      }
      else {
        theta.max <- theta
      }
      theta <- runif(n = 1,
                     min = theta.min,
                     max = theta.max)
    }
  }
  return(f.new)
}

