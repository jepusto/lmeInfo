
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lmeInfo

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/jepusto/lmeInfo.svg?branch=master)](https://travis-ci.org/jepusto/lmeInfo)
[![Codecov test
coverage](https://codecov.io/gh/jepusto/lmeInfo/branch/master/graph/badge.svg)](https://codecov.io/gh/jepusto/lmeInfo?branch=master)
<!-- badges: end -->

## Information matrices for fitted nlme::lme models

`lmeInfo` provides analytic derivatives and information matrices for
fitted lmeStruct models from `nlme::lme()`. The package includes
functions for estimating the sampling variance-covariance of variance
component parameters using the inverse Fisher information for fitted LME
models. The variance component paramters include the parameters of the
random effects structure, the variance structure and the correlation
structure. The package entails the calculation of the inverse of the
expected or average Fisher information matrix of the full maximum or
restricted maximum likelihood. Function is also available for
calculating design-comparable standardized mean difference effect size
(Pustejovsky et al., 2014) for multiple baseline across behaviors
designs, multi-level multiple baseline designs, alternating treatment
designs and meta-analyses.

## Installation

You can install the released version of lmeInfo from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("lmeInfo")
```

You can also install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jepusto/lmeInfo")
```

## Demonstration

We use the dataset from Bryant et al. (2016) to demonstrate how
`lmeInfo` estimates the sampling variance-covariance of variance
component parameters of fitted LME models and calculates the
design-comparable standardized mean difference effect size for a
three-level multiple baseline design.

### Sampling variance-covariance of variance component parameters

``` r
library(lmeInfo)
library(nlme)
data(Bryant2016)

Bryant2016_RML <- lme(fixed = outcome ~ treatment,
                      random = ~ 1 | school/case,
                      correlation = corAR1(0, ~ session | school/case),
                      data = Bryant2016)
```

The sampling variance-covariance of variance component parameters of
model `Bryant2016_RML` can be estimated using `varcomp_vcov()` in
`lmeInfo`. The sampling variance-covariance of variance component
parameters is calculated using the inverse expected Fisher information
when the type argument in `varcomp_vcov()` is specified as **expected**.
The inverse of the average Fisher information matrix is estimated when
the type argument is specified as **average**.

``` r
varcomp_vcov(Bryant2016_RML, type = "expected")
#>                             Tau.school.var((Intercept)) Tau.case.var((Intercept))    cor_params      sigma_sq
#> Tau.school.var((Intercept))                8.133900e+04             -2.140289e+04 -0.0100434219 -4.701744e+02
#> Tau.case.var((Intercept))                 -2.140289e+04              3.996251e+05 -8.2280433012 -3.815479e+05
#> cor_params                                -1.004342e-02             -8.228043e+00  0.0002154736  9.872975e+00
#> sigma_sq                                  -4.701744e+02             -3.815479e+05  9.8729751627  4.602107e+05
```

### The design-comparable standardized mean difference effect size

The design-comparable standardized mean difference effect size is
calculated using `g_mlm()` in `lmeInfo`. The standard errors of the
unadjusted or adjusted effect size are estimated with the inverse
expected Fisher information matrix when the infotype argument in
`g_mlm()` is specified as **expected**. The inverse average Fisher
information matrix is used when the infotype argument is specified as
**average**. The estimation method is built on the method described in
Pustejovsky, Hedges, & Shadish (2014).

``` r
Bryant2016_g <- g_mlm(Bryant2016_RML, p_const = c(0,1), r_const = c(1,1,0,1), infotype = "expected")
```

``` r
print(Bryant2016_g)
#>                           est    se
#> unadjusted effect size  0.481 0.122
#> adjusted effect size    0.463 0.118
#> degree of freedom      20.169
```

## References

Bryant, B. R., Bryant, D. P., Porterfield, J., Dennis, M. S., Falcomata,
T., Valentine, C., … & Bell, K. (2016). The effects of a Tier 3
intervention on the mathematics performance of second grade students
with severe mathematics difficulties. Journal of learning disabilities,
49(2), 176-188. <doi:10.1177/0022219414538516>

Pustejovsky, J. E., Hedges, L. V., & Shadish, W. R. (2014).
Design-comparable effect sizes in multiple baseline designs: A general
modeling framework. Journal of Educational and Behavioral Statistics, 39
(5), 368–393. <doi:10.3102/1076998614547577>
