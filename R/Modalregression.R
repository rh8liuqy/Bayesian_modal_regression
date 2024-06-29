#' Bayesian Modal Regression
#'
#' @export
#' @param formula a formula.
#' @param data a dataframe.
#' @param model a description of the error distribution. Can be one of "FG", "DTP" and "TPSC".
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return A `draw` object from the \strong{posterior} package.
#' @details
#' The Bayesian modal regression model based on the FG, DTP or TPSC distribution is defined as:
#' \deqn{Y_{i} = \mathbf{X}_{i} \boldsymbol{\beta} + e_{i},}
#' where \eqn{e_{i}} follows the FG, DTP or TPSC distribution.
#'
#' More details of the Bayesian modal regression model can be found at at Liu, Huang, and Bai (2024) <https://arxiv.org/pdf/2211.10776>.
#'
#' @example /inst/examples/Modalregression_example.R
#' @references
#' \insertRef{liu2022bayesian}{GUD}
modal_regression <- function(formula, data, model,...) {
  if (! inherits(formula, "formula")) {
    stop("The first argument must be a formula.")
  }
  match.arg(model,c("FG","DTP","TPSC"))
  # ensure data is a dataframe
  data <- as.data.frame(data)
  X <- model.matrix(object = formula, data = data)
  coltxt <- colnames(X)
  N <- nrow(X)
  P <- ncol(X)
  y <- as.numeric(data[,as.character(formula[[2]])])
  dat <- list(N = N,
              P = P,
              y = y,
              X = X)
  if (model == "TPSC") {
    sampling_output <- rstan::sampling(stanmodels$TPSC, data = dat,...)
  } else if (model == "DTP") {
    sampling_output <- rstan::sampling(stanmodels$DTP, data = dat,...)
  } else if(model == "FG") {
    sampling_output <- rstan::sampling(stanmodels$FG, data = dat,...)
  }
  draws_output <- posterior::as_draws(sampling_output)
  if (model == "TPSC") {
    posterior::variables(draws_output)[4:(4+length(coltxt)-1)] <- coltxt
  } else if (model == "DTP") {
    posterior::variables(draws_output)[5:(5+length(coltxt)-1)] <- coltxt
  } else if (model == "FG") {
    posterior::variables(draws_output)[4:(4+length(coltxt)-1)] <- coltxt
  }

  output <- draws_output

  return(output)
}
