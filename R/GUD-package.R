#' The 'GUD' package.
#'
#' @description This R package encompasses the probability density functions of three key distributions: the flexible Gumbel distribution, the double two-piece Student-t distribution, and the two-piece scale Student-t distribution, all belonging to the general unimodal distribution family, along with their corresponding sampling algorithms. Additionally, the package offers a function for Bayesian linear modal regression, leveraging these three distributions for model fitting.
#'
#' @docType package
#' @name GUD
"_PACKAGE"
#' @useDynLib GUD, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats model.matrix rbinom rgamma rnorm runif
#' @importFrom Rdpack reprompt
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package version 2.32.3. https://mc-stan.org
#' \insertRef{liu2022bayesian}{GUD}
NULL
