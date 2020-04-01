
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
functions for estimating the inverse of the expected or average Fisher
information matrix of the full maximum or restricted maximum likelihood.
It entails estimating the sampling variance of variance parameters,
including the parameters of the random effects structure, variance
structure and correlation structure. Function is also available for
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
`lmeInfo` estimates the inverse Fisher information matrices for fitted
LME models and calculates the design-comparable standardized mean
difference effect size for a three-level multiple baseline design.

  - The inverse of Fisher information matrices

<!-- end list -->

``` r
library(lmeInfo)
library(nlme)
data(Bryant2016)

Bryant2016_RML <- lme(fixed = outcome ~ treatment,
                      random = ~ 1 | school/case,
                      correlation = corAR1(0, ~ session | school/case),
                      data = Bryant2016)
```

The inverse of the expected Fisher information matrix of model
`Bryant2016_RML` can be estimated using `Fisher_info()` in `lmeInfo`
when the type argument is specified as **expected**. The inverse of the
average Fisher information matrix is estimated when type is specified as
**average**.

``` r
Fisher_info(Bryant2016_RML, type = "expected")
#>                             Tau.school.var((Intercept)) Tau.case.var((Intercept))    cor_params      sigma_sq
#> Tau.school.var((Intercept))                1.322700e-05              3.481279e-06  4.031958e-02  2.034762e-06
#> Tau.case.var((Intercept))                  3.481279e-06              1.299374e-05  1.510332e-01  7.536162e-06
#> cor_params                                 4.031958e-02              1.510332e-01  2.744559e+05 -5.762687e+00
#> sigma_sq                                   2.034762e-06              7.536162e-06 -5.762687e+00  1.320509e-04
```

  - The design-comparable standardized mean difference effect size

<!-- end list -->

``` r
Bryant2016_g <- g_mlm(Bryant2016_RML, p_const = c(0,1), r_const = c(1,1,0,1))
```

``` r
print(Bryant2016_g)
#>                           est    se
#> unadjusted effect size  0.481 0.122
#> adjusted effect size    0.463 0.118
#> degree of freedom      20.169
```

Pustejovsky, Hedges & Shadish (2014) provides detailed estimation
methods for the effect sizes, standard errors and degree of freedom.

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
