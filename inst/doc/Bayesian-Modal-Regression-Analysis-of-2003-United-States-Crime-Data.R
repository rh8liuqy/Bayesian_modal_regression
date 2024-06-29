## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(GUD)

## ----fig.width=6,fig.height=6-------------------------------------------------
# load data crime from the GUD package
df1 <- crime
# the conditional scatter plot matrices of U.S. crime data
if (require(lattice)) {
  lattice::splom(~df1[c(6,4,9,3)],
                 main = NULL,
                 panel = function(x,y,...) {
                   panel.splom(x,y,...)
            })
}

## -----------------------------------------------------------------------------
TPSC_model <- modal_regression(`murder rate` ~ college + poverty + metropolitan, 
                               data = df1, 
                               model = "TPSC",
                               chains = 2,
                               iter = 2000)

## -----------------------------------------------------------------------------
print(summary(TPSC_model), n = 7)

## ----fig.width=6, fig.height=4------------------------------------------------
if (require(bayesplot)) {
  bayesplot::mcmc_trace(TPSC_model, pars = c("(Intercept)",
                                             "college", 
                                             "poverty", 
                                             "metropolitan"))
}

## -----------------------------------------------------------------------------
summary(posterior::subset_draws(TPSC_model, variable = "ystar"))

