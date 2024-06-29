#' The Density Function of the Gumbel Distribution
#'
#' @param x vector of quantiles.
#' @param loc vector of location parameters.
#' @param sigma vector of scale parameters.
#'
#' @return vector of densities.
#' @noRd
dgumbel <- function(x,loc,sigma){
  z <- (x-loc)/sigma
  return(1/sigma*exp(-z-exp(-z)))
}

#' The Random Number Generation of the Gumbel Distribution
#'
#' @param n number of observations.
#' @param loc vector of location parameters.
#' @param sigma vector of scale parameters.
#' @return vector of random deviates.
#' @noRd
rgumbel <- function(n,loc,sigma){
  U <- runif(n = n)
  output <- loc - sigma * log(-log(U))
  return(output)
}

#' The Flexible Gumbel Distribution
#' @param x vector of quantiles.
#' @param n number of observations.
#' @param w vector of weight parameters.
#' @param loc vector of the location parameters.
#' @param sigma1 vector of the scale parameters of the left skewed part.
#' @param sigma2 vector of the scale parameters of the right skewed part.
#' @name dFG
#' @details
#' The Gumbel distribution has the density
#' \deqn{f_{\text {Gumbel }}(y \mid \theta, \sigma)=\frac{1}{\sigma} \exp \left\{-\frac{y-\theta}{\sigma}-\exp \left(-\frac{y-\theta}{\sigma}\right)\right\},}
#' where \eqn{\theta \in \mathbb{R}} is the mode as the location parameter, \eqn{\sigma > 0} is the scale parameter.
#'
#' The flexible Gumbel distribution has the density
#' \deqn{f_{\mathrm{FG}}\left(y \mid w, \theta, \sigma_1, \sigma_2\right)=w f_{\text {Gumbel }}\left(-y \mid-\theta, \sigma_1\right)+(1-w) f_{\text {Gumbel }}\left(y \mid \theta, \sigma_2\right) .}
#' where \eqn{w \in [0,1]} is the weight parameter, \eqn{\sigma_{1} > 0} is the scale parameter of the left skewed part and \eqn{\sigma_{2} > 0} is the scale parameter of the right skewed part.
#' @return `dFG` gives the density. `rFG` generates random deviates.
#' @references
#' \insertRef{liu2022bayesian}{GUD}
#' @examples
#' set.seed(100)
#' require(graphics)
#'
#' # Random Number Generation
#' X <- rFG(n = 1e5, w = 0.3, loc = 0, sigma1 = 1, sigma2 = 2)
#'
#' # Plot the histogram
#' hist(X, breaks = 100, freq = FALSE)
#'
#' # The red dashed line should match the underlining histogram
#' points(x = seq(-10,20,length.out = 1000),
#'        y = dFG(x = seq(-10,20,length.out = 1000),
#'                w = 0.3, loc = 0, sigma1 = 1, sigma2 = 2),
#'        type = "l",
#'        col = "red",
#'        lwd = 3,
#'        lty = 2)
NULL
#> NULL

#' @rdname dFG
#' @export
dFG <- function(x,w,loc,sigma1,sigma2) {
  return(return(w*dgumbel(-x,-loc,sigma1)+(1-w)*dgumbel(x,loc,sigma2)))
}

#' @rdname dFG
#' @export
rFG <- function(n,w,loc,sigma1,sigma2) {
  Z <- rbinom(n = n, size = 1, prob = w)
  X1 <- - rgumbel(n = n,loc = 0,sigma = sigma1)
  X2 <- rgumbel(n = n,loc = 0,sigma = sigma2)
  output <- Z * X1 + (1 - Z) * X2
  return(output)
}
