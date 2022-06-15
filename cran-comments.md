## Resubmit comments

In this version we have:
* Added an option to allow heterogeneous level-1 variance components.
* Generalized the `g_mlm()` function to allow use of separate models for the numerator and denominator of the effect size.
* Modified the stored results of `g_mlm` so that the `returnModel` argument is no longer necessary.
* Fixed a bug in handling of models with missing observations.

## Test environments
* local Windows 10 Pro, R 4.1.0
- ubuntu 20.04.4 LTS (on Github), R devel, release, oldrelease
- macOS-latest (on Github), R release
- Windows-latest (on Github), R release
* win-builder (devel, release, oldrelease)
* r-hub:
  - Debian Linux, R-devel, clang, ISO-8859-15 locale
  - Debian Linux, R-devel, GCC
  - Debian Linux, R-patched, GCC
  - Debian Linux, R-release, GCC
  - Fedora Linux, R-devel, clang, gfortran
  - Fedora Linux, R-devel, GCC
  - macOS 10.13.6 High Sierra, R-release, CRAN's setup
  - Apple Silicon (M1), macOS 11.6 Big Sur, R-release
  - Windows Server 2022, R-devel, 64 bit
  - Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit
  - Windows Server 2008 R2 SP1, R-release, 32/64 bit

## R CMD check results
There was 1 NOTE:

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1177/0022219414538516
    From: man/Bryant2016.Rd
    Status: 503
    Message: Service Unavailable
  URL: https://doi.org/10.2307/2533274
    From: inst/doc/Information-matrices-for-fitted-LME-models.html
    Status: 403
    Message: Forbidden
  URL: https://doi.org/10.2307/2533558
    From: inst/doc/Information-matrices-for-fitted-LME-models.html
    Status: 403
    Message: Forbidden
  URL: https://doi.org/10.1080/01621459.1988.10478693
    From: inst/doc/Information-matrices-for-fitted-LME-models.html
    Status: 403
    Message: Forbidden
  URL: https://doi.org/10.3102/1076998614547577
    From: man/g_mlm.Rd
    Status: 503
    Message: Service Unavailable

Found the following (possibly) invalid DOIs:
  DOI: 10.3102/1076998614547577
    From: DESCRIPTION
    Status: Service Unavailable
    Message: 503

## revdepcheck results
We checked 3 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

* We saw 0 new problems
* We failed to check 0 packages
