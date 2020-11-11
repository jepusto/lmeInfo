## Resubmission
This is a resubmission. This resubmission corrects the errors shown on CRAN package checks. We fix a bug leading to errors in `Fisher_info()` with models that have more than two levels. We add tests for handling models with missing outcome and/or covariate values.

## Test environments
* local Windows 10 Pro, R 4.0.2
* ubuntu 16.04.6 LTS (on travis-ci), R-release, R-devel
* win-builder (devel, release, oldrelease)
* r-hub:
  * Debian Linux, R-devel, clang, ISO-8859-15 locale
  * Debian Linux, R-devel, GCC
  * Debian Linux, R-patched, GCC
  * Debian Linux, R-release, GCC
  * Fedora Linux, R-devel, clang, gfortran
  * Fedora Linux, R-devel, GCC
  * macOS 10.13.6 High Sierra, R-release, CRAN's setup
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  * Windows Server 2008 R2 SP1, R-oldrel, 32/64 bit
  * Windows Server 2008 R2 SP1, R-release, 32/64 bit

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.

