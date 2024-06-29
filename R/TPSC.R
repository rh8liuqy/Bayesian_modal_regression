#' The TPSC-Student-t Distribution
#' @param x vector of quantiles.
#' @param n number of observations.
#' @param w vector of weight parameters.
#' @param theta vector of the location parameters.
#' @param sigma vector of the scale parameters.
#' @param delta the degree of freedom.
#' @name dTPSC
#' @details
#' The TPSC-Student-t distribution has the density
#' \deqn{f_{\mathrm{TPSC}}(y \mid w, \theta, \sigma, \delta)=w f_{\mathrm{LT}}\left(y \mid \theta, \sigma \sqrt{\frac{w}{1-w}}, \delta\right)+(1-w) f_{\mathrm{RT}}\left(y \mid \theta, \sigma \sqrt{\frac{1-w}{w}}, \delta\right),}
#' where \deqn{f_{\mathrm{LT}}(y \mid \theta, \sigma, \delta)=\frac{2}{\sigma} f\left(\left.\frac{y-\theta}{\sigma} \right\rvert\, \delta\right) \mathbb{I}(y<\theta),} and \deqn{f_{\mathrm{RT}}(y \mid \theta, \sigma, \delta)=\frac{2}{\sigma} f\left(\left.\frac{y-\theta}{\sigma} \right\rvert\, \delta\right) \mathbb{I}(y \geq \theta).}
#' Additionally, \eqn{f(y \mid \delta)} represents the density function of the standardized Student-t distribution with the degree of freedom \eqn{\delta}.
#'
#' @return `dTPSC` gives the density. `rTPSC` generates random deviates.
#' @references
#' \insertRef{liu2022bayesian}{GUD}
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
