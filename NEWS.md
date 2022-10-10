# lmeInfo 0.3.0

* Fixed a bug in the `Fisher_info()` function that causes incorrect order for the three-level models and for the models with multiple random effects at the same level.
* Removed `cnvg_warn` argument from the `g_mlm()` function.
* Updated the unit tests related to Lambert et al. (2006) data.

# lmeInfo 0.2.1

* Fixed a bug in `extract_varcomp()` that caused some variance components to be dropped if the variables involved in the random effects formula involved special characters such as `.`, `(`, `)`, or `^`.

# lmeInfo 0.2.0

* Added an option to return separate level-1 variance components for models that use `weights = varIdent(form = ~ 1 | Stratum)`.
* Generalized the `g_mlm()` function to allow use of separate models for the numerator and denominator of the effect size.
* Modified the stored results of `g_mlm()` so that the `returnModel` argument is no longer necessary.
* Fixed a bug in handling of models with missing observations that have `na.action` of `na.exclude()`.
* Fixed a bug in internal functions for constructing level-1 variance covariance structures in models with multi-variate structure (i.e., models with non-null `modelStruct$corStruct`).

# lmeInfo 0.1.2

* Corrected a bug leading to errors in `Fisher_info()` with models that have more than two levels.
* Added informative errors for `summary.g_mlm()`. `summary()` method only works for a `g_mlm()` object when setting `returnModel = TRUE`, which is the default.
* Updated the formula for standard error of effect size in the vignette to match the method in `g_mlm()`.
* Added tests for handling models with missing outcome and/or covariate values.

# lmeInfo 0.1.1

* Fixed the "Additional issues" in the unit tests identified by the CRAN package checks. 
* Added package website.

# lmeInfo 0.1.0

* Initial release.
