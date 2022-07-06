## Resubmit comments

In this version we have fixed a bug in `extract_varcomp()` that caused some variance components to be dropped if the variables involved in the random effects formula involved special characters such as `.`, `(`, `)`, or `^`. This bug is important to correct for the package's downstream dependency (the scdhlm package). We apologize for the inconvenience caused by our rapid re-submission. Thank you.

## Test environments

* local Windows 10 Pro, R 4.1.0
* ubuntu 20.04.4 LTS (on Github), R devel, release, oldrelease
* macOS-latest (on Github), R release
* Windows-latest (on Github), R release
* win-builder (devel, release, oldrelease)
* r-hub:
  - macOS 10.13.6 High Sierra, R-release, CRAN's setup
  - Apple Silicon (M1), macOS 11.6 Big Sur, R-release

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
