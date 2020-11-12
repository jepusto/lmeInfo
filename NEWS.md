# lmeInfo 0.1.2

* Fixed the errors in the unit tests of `check_against_scdhlm()` detected by the CRAN package checks.
* Corrected a bug leading to errors in `Fisher_info()` with models that have more than two levels.
* Added tests for handling models with missing outcome and/or covariate values.
* Added informative errors for `summary.g_mlm()`. `summary()` method only works for a `g_mlm()` object when setting `returnModel = TRUE`, which is the default.
* Updated the formula for standard error of effect size in the vignette to match the method in `g_mlm()`.

# lmeInfo 0.1.1

* Fixed the "Additional issues" in the unit tests identified by the CRAN package checks. 
* Added package website.

# lmeInfo 0.1.0

* Initial release.
