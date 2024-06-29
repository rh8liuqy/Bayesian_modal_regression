
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Installation

The `GUD`
[![](https://cranlogs.r-pkg.org/badges/grand-total/GUD)](https://cran.r-project.org/web/packages/GUD/index.html)
package provides probability density functions and sampling algorithms
for three key distributions from the General Unimodal Distribution (GUD)
family: the Flexible Gumbel (FG) distribution, the Double Two-Piece
(DTP) Student-t distribution, and the Two-Piece Scale (TPSC) Student-t
distribution. Additionally, this package includes a function for
Bayesian linear modal regression, leveraging these three distributions
for model fitting.

You can install the `GUD` package:

``` r
if (! require(GUD)) {
  install.packages(pkgs = "GUD",
                   repos = "https://cran.rstudio.com/" )
} else {
  require(GUD)
}
#> Loading required package: GUD
```

## Introduction

The general unimodal distribution (GUD) family is essentially a family
of two-component mixture distributions. The probability density function
(pdf) of a member of GUD family is $$
f\left(y \mid w, \theta, \boldsymbol{\xi}_1, \boldsymbol{\xi}_2\right)=w f_1\left(y \mid \theta, \boldsymbol{\xi}_1\right)+(1-w) f_2\left(y \mid \theta, \boldsymbol{\xi}_2\right),
$$ where $w \in [0,1]$ is the weight parameter,
$\theta \in (-\infty, +\infty)$ is the mode as a location parameter,
$\boldsymbol{\xi}_1$ consists of parameters other than the location
parameter in $f_1\left(\cdot \mid \theta, \boldsymbol{\xi}_1\right)$ and
$\boldsymbol{\xi}_2$ is defined similarly for
$f_2\left(\cdot \mid \theta, \boldsymbol{\xi}_2\right)$. Besides
unimodality, all members of the GUD family share three features:

1.  The pdfs $f_1\left(\cdot \mid \theta, \boldsymbol{\xi}_1\right)$ and
    $f_2\left(\cdot \mid \theta, \boldsymbol{\xi}_2\right)$ are unimodal
    at $\theta$.
2.  The pdfs $f_1\left(\cdot \mid \theta, \boldsymbol{\xi}_1\right)$ and
    $f_2\left(\cdot \mid \theta, \boldsymbol{\xi}_2\right)$ are
    left-skewed and right-skewed respectively.
3.  The mixture pdf
    $f\left(\cdot \mid w, \theta, \boldsymbol{\xi}_1, \boldsymbol{\xi}_2\right)$
    in (1) is continuous in its domain.

More details of the GUD family can be found in [Liu, Q., Huang, X., &
Bai, R. (2024)](https://arxiv.org/abs/2211.10776).

## Bayesian Modal Regression Analysis of 2003 United States Crime Data

In this section, we demonstrate how to use the `GUD` package to analyze
2003 United States Crime Data as in Section 2 of [Liu, Q., Huang, X., &
Bai, R. (2024)](https://arxiv.org/abs/2211.10776).

In “The Art and Science of Learning from Data, 5th edition” by Alan
Agresti, Christine A. Franklin, and Bernhard Klingenberg, an interesting
example about the 2003 United States crime data is presented to
demonstrate the influence of outliers in the classic linear regression
model. This example is very compelling and partially motivates the
construction of the Bayesian modal regression based on the GUD family.
This data contains the murder rate, percentage of college education,
poverty percentage, and metropolitan rate for the 50 states in the
United States and the District of Columbia (D.C.) from 2003. The murder
rate is defined as the annual number of murders per $100{,}000$ people
in the population. The poverty percentage is the percentage of residents
with income below the poverty level, and the metropolitan rate is
defined as the percentage of population living in the metropolitan area.
In the exploratory data analysis, we present the conditional scatter
plot matrices below.

``` r
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
#> Loading required package: lattice
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

In the conditional scatter plot matrices, we notice an outlier,
Washington, D.C., which stands out and does not follow the common
pattern of other states.

Next, we demonstrate how to fit the Bayesian modal regression model
based on the TPSC distribution to the 2003 United States crime data.

``` r
TPSC_model <- modal_regression(`murder rate` ~ college + poverty + metropolitan, 
                               data = df1, 
                               model = "TPSC",
                               chains = 2,
                               iter = 2000)
#> 
#> SAMPLING FOR MODEL 'TPSC' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3.1e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.31 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 1.128 seconds (Warm-up)
#> Chain 1:                0.734 seconds (Sampling)
#> Chain 1:                1.862 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'TPSC' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2.1e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.21 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 1.054 seconds (Warm-up)
#> Chain 2:                0.799 seconds (Sampling)
#> Chain 2:                1.853 seconds (Total)
#> Chain 2:
```

### Summary of Bayesian Analysis

One can summarize the Bayesian analysis using the summary function.

``` r
print(summary(TPSC_model), n = 7)
#> # A tibble: 113 × 10
#>   variable        mean  median     sd    mad       q5     q95  rhat ess_bulk
#>   <chr>          <dbl>   <dbl>  <dbl>  <dbl>    <dbl>   <dbl> <dbl>    <dbl>
#> 1 w             0.270   0.276  0.119  0.134   0.0790   0.466   1.01     347.
#> 2 delta         1.89    1.80   0.645  0.562   1.08     3.01    1.00    1024.
#> 3 sigma         1.17    1.16   0.259  0.254   0.769    1.62    1.00     464.
#> 4 (Intercept)   1.11    1.13   2.65   2.69   -3.53     5.42    1.00     706.
#> 5 college      -0.199  -0.201  0.0804 0.0824 -0.328   -0.0622  1.00     770.
#> 6 poverty       0.240   0.248  0.135  0.131   0.00362  0.446   1.01     503.
#> 7 metropolitan  0.0640  0.0623 0.0158 0.0156  0.0397   0.0916  1.01     489.
#> # ℹ 106 more rows
#> # ℹ 1 more variable: ess_tail <dbl>
```

One can present the traceplot of the MCMC chain using the
`bayesplot::mcmc_trace` function.

``` r
if (require(bayesplot)) {
  bayesplot::mcmc_trace(TPSC_model, pars = c("(Intercept)",
                                             "college", 
                                             "poverty", 
                                             "metropolitan"))
}
#> Loading required package: bayesplot
#> Warning: package 'bayesplot' was built under R version 4.3.1
#> This is bayesplot version 1.11.1
#> - Online documentation and vignettes at mc-stan.org/bayesplot
#> - bayesplot theme set to bayesplot::theme_default()
#>    * Does _not_ affect other ggplot2 plots
#>    * See ?bayesplot_theme_set for details on theme setting
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

The summary of posterior predictive distribution can be assessed using
the following command. Here `ystar[1]` represents the posterior
prediction of the first observation in the dataset, and so on.

``` r
summary(posterior::subset_draws(TPSC_model, variable = "ystar"))
#> # A tibble: 51 × 10
#>    variable   mean median    sd   mad     q5   q95  rhat ess_bulk ess_tail
#>    <chr>     <dbl>  <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl>    <dbl>    <dbl>
#>  1 ystar[1]   7.86   6.12 20.0   2.04  3.05  13.8  0.999    1974.    1853.
#>  2 ystar[2]   2.31   1.28  6.66  2.02 -1.69   8.32 1.00     1860.    1872.
#>  3 ystar[3]   7.16   6.12 15.5   1.96  3.33  13.9  0.999    1933.    1897.
#>  4 ystar[4]   6.56   5.61 16.8   2.21  2.27  11.9  1.00     1680.    1791.
#>  5 ystar[5]   6.87   6.14  6.12  1.94  3.17  13.0  0.999    1969.    1961.
#>  6 ystar[6]   4.14   2.71 13.8   2.05 -0.411 11.0  1.00     1641.    1829.
#>  7 ystar[7]   5.21   3.74 16.8   1.95  0.746 11.6  1.00     1981.    2096.
#>  8 ystar[8]   5.79   4.87 12.0   1.89  2.24  11.8  1.00     1808.    1924.
#>  9 ystar[9]   7.24   5.21 49.4   2.74  1.56  13.6  1.00     1421.    1758.
#> 10 ystar[10]  8.15   6.42 38.9   2.00  3.48  13.3  1.00     2092.    1963.
#> # ℹ 41 more rows
```

Further comparisons between mean, median, and modal regression can be
found in Section 2 of [Liu, Q., Huang, X., & Bai, R.
(2024)](https://arxiv.org/abs/2211.10776) and Section 6 of [Liu, Q.,
Huang, X., & Zhou, H. (2024)](https://www.mdpi.com/2571-905X/7/1/19).
