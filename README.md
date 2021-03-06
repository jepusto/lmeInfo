
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lmeInfo

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/jepusto/lmeInfo.svg?branch=master)](https://travis-ci.org/jepusto/lmeInfo)
[![Codecov test
coverage](https://codecov.io/gh/jepusto/lmeInfo/branch/master/graph/badge.svg)](https://codecov.io/gh/jepusto/lmeInfo?branch=master)
[![](http://www.r-pkg.org/badges/version/lmeInfo)](https://CRAN.R-project.org/package=lmeInfo)
[![](http://cranlogs.r-pkg.org/badges/grand-total/lmeInfo)](https://CRAN.R-project.org/package=lmeInfo)
[![](http://cranlogs.r-pkg.org/badges/last-month/lmeInfo)](https://CRAN.R-project.org/package=lmeInfo)
<!-- badges: end -->

## Information matrices for fitted `lme` and `gls` models

`lmeInfo` provides analytic derivatives and information matrices for
fitted linear mixed effects (lme) models and generalized least squares
(gls) models estimated using `nlme::lme()` and `nlme::gls()`,
respectively. The package includes functions for estimating the sampling
variance-covariance of variance component parameters using the inverse
Fisher information. The variance components include the parameters of
the random effects structure (for lme models), the variance structure,
and the correlation structure. The expected and average forms of the
Fisher information matrix are used in the calculations, and models
estimated by full maximum likelihood or restricted maximum likelihood
are supported. The package also includes a function for estimating
standardized mean difference effect sizes (Pustejovsky et al., 2014)
based on fitted `lme` or `gls` models.

## Installation

You can install the released version of lmeInfo from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("lmeInfo")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jepusto/lmeInfo")
```

## Demonstration

We use a dataset from a multiple baseline study conducted by Bryant and
colleagues (2016) to demonstrate how `lmeInfo` estimates the sampling
variance-covariance of variance component parameters of fitted LME
models. We also show how to calculate a design-comparable standardized
mean difference effect size based on the fitted model.

The study by Bryant and colleagues (2016) involved collecting repeated
measures of math performance on multiple students, in each of three
schools. After an initial baseline period in each school, an
intervention was introduced and its effects on student math performance
were observed over time. For sake of illustration, we use a very simple
model for these data, consisting of a simple change in levels coinciding
with the introduction of treatment. We include random effects for each
school and each student. Here we fit the model using `nlme::lme()`:

``` r
library(lmeInfo)
library(nlme)
data(Bryant2016)

Bryant2016_RML <- lme(fixed = outcome ~ treatment,
                      random = ~ 1 | school/case,
                      data = Bryant2016)

summary(Bryant2016_RML)
#> Linear mixed-effects model fit by REML
#>  Data: Bryant2016 
#>        AIC      BIC    logLik
#>   2627.572 2646.041 -1308.786
#> 
#> Random effects:
#>  Formula: ~1 | school
#>         (Intercept)
#> StdDev:    12.57935
#> 
#>  Formula: ~1 | case %in% school
#>         (Intercept) Residual
#> StdDev:    15.98217   18.398
#> 
#> Fixed effects: outcome ~ treatment 
#>                       Value Std.Error  DF   t-value p-value
#> (Intercept)        56.14333  9.019268 286  6.224821       0
#> treatmenttreatment 49.33454  2.399065 286 20.564070       0
#>  Correlation: 
#>                    (Intr)
#> treatmenttreatment -0.194
#> 
#> Standardized Within-Group Residuals:
#>         Min          Q1         Med          Q3         Max 
#> -3.15430670 -0.65412204  0.02112957  0.61701090  2.90624593 
#> 
#> Number of Observations: 299
#> Number of Groups: 
#>           school case %in% school 
#>                3               12
```

The estimated variance components from the fitted model can be obtained
using `extract_varcomp()`:

``` r
extract_varcomp(Bryant2016_RML)
#> $Tau
#> $Tau$school
#> school.var((Intercept)) 
#>                  158.24 
#> 
#> $Tau$case
#> case.var((Intercept)) 
#>              255.4297 
#> 
#> 
#> $cor_params
#> numeric(0)
#> 
#> $var_params
#> numeric(0)
#> 
#> $sigma_sq
#> [1] 338.4864
#> 
#> attr(,"class")
#> [1] "varcomp"
```

### Sampling variance-covariance of variance component parameters

The sampling variance-covariance of variance component parameters of
model `Bryant2016_RML` can be estimated using `varcomp_vcov()` in
`lmeInfo`. Setting `type = "expected"` will calculate the sampling
variance-covariance of variance component parameters using the inverse
expected Fisher information. Setting `type = "average"` will calculate
the inverse of the average Fisher information matrix (Gilmour, Thompson,
& Cullis, 1995).

``` r
varcomp_vcov(Bryant2016_RML, type = "expected")
#>                             Tau.school.var((Intercept)) Tau.case.var((Intercept))     sigma_sq
#> Tau.school.var((Intercept))                5.695541e+04               -4547.41772   0.06132317
#> Tau.case.var((Intercept))                 -4.547418e+03               16045.57184 -32.27967543
#> sigma_sq                                   6.132317e-02                 -32.27968 801.20973947
```

### Estimating a standardized mean difference effect size

The package also includes a function, `g_mlm()`, for estimating a
standardized mean difference effect size from a fitted multi-level
model. The estimation methods follow Pustejovsky, Hedges, and Shadish
(2014). A standardized mean difference effect size parameter can be
defined as the ratio of a linear combination of the model’s fixed effect
parameters over the square root of a linear combination of the model’s
variance components. The `g_mlm()` function takes as inputs a fitted
multi-level model and the vectors `p_const` and `r_const`, which define
the linear combinations of fixed effects and variance components,
respectively. The function calculates an effect size estimate by first
substituting maximum likelihood or restricted maximum likelihood
estimates in place of the corresponding parameters, then applying a
small-sample correction. The small-sample correction and the standard
error are based on approximating the distribution of the estimator by a
t distribution, with degrees of freedom given by a Satterthwaite
approximation (Pustejovsky, Hedges, & Shadish, 2014). The `g_mlm()`
function includes an option allowing use of the expected or average form
of the Fisher information matrix in the calculations. The `g_mlm()`
function also includes an option allowing returning the fitted model
parameters in addition to effect size estimate.

In our model for the Bryant data, we use the treatment effect in the
numerator of the effect size and the sum of the school-level,
student-level, and within-student variance components in the denominator
of the effect size. The constants are therefore given by `p_const =
c(0, 1)` and `r_const = c(1, 1, 1)`. The effect size estimate can be
calculated as:

``` r
Bryant2016_g <- g_mlm(Bryant2016_RML, p_const = c(0,1), r_const = c(1,1,1), infotype = "expected")
```

``` r
Bryant2016_g
#>                           est    se
#> unadjusted effect size  1.799 0.340
#> adjusted effect size    1.721 0.325
#> degree of freedom      17.504
```

A `summary()` method is also included, which includes more detail about
the model parameter estimates and effect size estimate when setting
`returnModel = TRUE` (the default) in `g_mlm()`:

``` r
summary(Bryant2016_g)
#>                                            est      se
#> Tau.school.school.var((Intercept))     158.240 238.653
#> Tau.case.case.var((Intercept))         255.430 126.671
#> sigma_sq                               338.486  28.306
#> total variance                         752.156 254.250
#> (Intercept)                             56.143   9.019
#> treatmenttreatment                      49.335   2.399
#> treatment effect at a specified time    49.335   2.399
#> unadjusted effect size                   1.799   0.340
#> adjusted effect size                     1.721   0.325
#> degree of freedom                       17.504        
#> constant kappa                           0.087        
#> logLik                               -1308.786
```

## Supported model components

The `Fisher_info()` and `varcomp_vcov()` functions operate on fitted
`lme()` and `gls()` models. Models fitted using `lme()` can include
three types of variance component parameters: random effects variances
and covariances, correlation structure parameters, and variance
structure parameters. Models fitted using `gls()` can include
correlation structure parameters and variance structure parameters. The
`nlme` package provides many different forms for each of these
components, not all of which are supported in `lmeInfo`. The package can
handle the following classes of variance components:

  - Random effects structure
      - `pdSymm` matrices, including in the `pdLogChol` and `pdNatural`
        parameterizations
      - `pdDiag` matrices
  - Correlation structure
      - `corAR1`
      - `corCAR1`
      - `corARMA` for MA(1) models only
      - `corCompSymm`
      - `corSymm`
  - Variance structure
      - `varIdent`
      - `varExp`
      - `varPower`
      - `varConstPower`

Calling `Fisher_info()` or `varcomp_vcov()` on a fitted model that
includes variance component structures outside of the supported classes
will trigger an informative error message.

## Related work

The `lmeInfo` package enhances the functionality of the [`nlme`
package](https://cran.r-project.org/package=nlme). However, it does not
work on `nlme()` models. The [`merDeriv`
package](https://CRAN.R-project.org/package=merDeriv) (Wang & Merkle,
2018) provides some related functionality for linear mixed effects
models estimated with `lme4::lmer()`, including variance-covariance
matrices for estimaed variance components based on the inverse of the
expected or observed Fisher information. Currently, the functionality in
`merDeriv` is limited to models with only one level of random effects.

## References

Bryant, B. R., Bryant, D. P., Porterfield, J., Dennis, M. S., Falcomata,
T., Valentine, C., … & Bell, K. (2016). The effects of a Tier 3
intervention on the mathematics performance of second grade students
with severe mathematics difficulties. *Journal of Learning Disabilities,
49*(2), 176-188. <https://doi.org/10.1177/0022219414538516>

Gilmour, A. R., Thompson, R., & Cullis, B. R. (1995). Average
information REML: An efficient algorithm for variance parameter
estimation in linear mixed models. *Biometrics, 51*(4), 1440–1450.
<https://doi.org/10.2307/2533274>

Pustejovsky, J. E., Hedges, L. V., & Shadish, W. R. (2014).
Design-comparable effect sizes in multiple baseline designs: A general
modeling framework. *Journal of Educational and Behavioral Statistics,
39*(5), 368–393. <https://doi.org/10.3102/1076998614547577>

Wang T, Merkle EC (2018). merDeriv: Derivative Computations for Linear
Mixed Effects Models with Application to Robust Standard Errors.
*Journal of Statistical Software, Code Snippets*, *87*(1), 1-16.
<https://doi.org/10.18637/jss.v087.c01>.
