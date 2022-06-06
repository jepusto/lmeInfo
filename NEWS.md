# lmeInfo 0.1.3.9999

* Added an option to return separate level-1 variance components for models that use `weights = varIdent(form = ~ 1 | Stratum)`.
* Generalized the `g_mlm()` function to allow use of separate models for the numerator and denominator of the effect size.
* Modified the stored results of `g_mlm()` so that the `returnModel` argument is no longer necessary.

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
