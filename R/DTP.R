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

#' The DTP-Student-t Distribution
#' @param x vector of quantiles.
#' @param n number of observations.
#' @param theta vector of the location parameters.
#' @param sigma1 vector of the scale parameters of the left skewed part.
#' @param sigma2 vector of the scale parameters of the right skewed part.
#' @param delta1 the degree of freedom of the left skewed part.
#' @param delta2 the degree of freedom of the right skewed part.
#' @name dDTP
#' @details
#' The DTP-Student-t distribution has the density
#' \deqn{f_{\mathrm{DTP}}\left(y \mid \theta, \sigma_1, \sigma_2, \delta_1, \delta_2\right)=w f_{\mathrm{LT}}\left(y \mid \theta, \sigma_1, \delta_1\right)+(1-w) f_{\mathrm{RT}}\left(y \mid \theta, \sigma_2, \delta_2\right),}
#' where \deqn{w=\frac{\sigma_1 f\left(0 \mid \delta_2\right)}{\sigma_1 f\left(0 \mid \delta_2\right)+\sigma_2 f\left(0 \mid \delta_1\right)},}
#' \eqn{f(0 \mid \delta)} represents \deqn{f((y-\theta) / \sigma \mid \delta)\text{ evaluated at } y=\theta,}
#' \deqn{f_{\mathrm{LT}}(y \mid \theta, \sigma, \delta)=\frac{2}{\sigma} f\left(\left.\frac{y-\theta}{\sigma} \right\rvert\, \delta\right) \mathbb{I}(y<\theta),} and
#' \deqn{f_{\mathrm{RT}}(y \mid \theta, \sigma, \delta)=\frac{2}{\sigma} f\left(\left.\frac{y-\theta}{\sigma} \right\rvert\, \delta\right) \mathbb{I}(y \geq \theta).}
#' Additionally, \eqn{f(y \mid \delta)} represents the density function of the standardized Student-t distribution with the degree of freedom \eqn{\delta}.
#'
#' @return `dDTP` gives the density. `rDTP` generates random deviates.
#' @references
#' \insertRef{liu2022bayesian}{GUD}
#' @example /inst/examples/DTP_example.R
NULL
#> NULL

#' @rdname dDTP
#' @export
dDTP <- function(x,theta,sigma1,sigma2,delta1,delta2){
  w <- sigma1*dtdist(0,delta2)/(sigma1*dtdist(0,delta2)+sigma2*dtdist(0,delta1))
  p1 <- w*2/sigma1*dtdist((x-theta)/sigma1,delta1)*(x<theta)
  p2 <- (1-w)*2/sigma2*dtdist((x-theta)/sigma2,delta2)*(x>=theta)
  return(p1+p2)
}

#' @rdname dDTP
#' @export
rDTP <- function(n,theta,sigma1,sigma2,delta1,delta2){
  w <- sigma1 * dtdist(x = 0,delta = delta2) / (sigma1 * dtdist(x = 0,delta = delta2) + sigma2 * dtdist(x = 0,delta = delta1))
  X1 <- -abs(rtdist(n = n,delta = delta1)) * sigma1
  X2 <- abs(rtdist(n = n,delta = delta2)) * sigma2
  Z <- rbinom(n = n, size = 1, prob = w)
  output <- X1*Z + X2*(1-Z) + theta
  return(output)
}
