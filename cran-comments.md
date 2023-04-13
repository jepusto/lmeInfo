## Resubmit comments

This is a resubmission. This resubmission fixes the additional issues (ATLAS) in CRAN package checks.

## Test environments

* local Windows 10 Pro, R 4.2.2
* ubuntu 20.04.5 LTS (on Github), R devel, release, oldrelease
* macOS-latest (on Github), R release
* Windows-latest (on Github), R release
* win-builder (devel, release, oldrelease)
* r-hub:
  * Windows Server 2022, R-devel, 64 bit
  * Windows Server 2022, R-oldrel, 32/64 bit
  * Fedora Linux, R-devel, clang, gfortran


## R CMD check results
There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Possibly mis-spelled words in DESCRIPTION:
    Pustejovsky (24:48)
    Shadish (24:73)
    gls (16:77, 25:28)
    lme (16:34, 20:57, 25:21)
  
  These words are spelled correctly.

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

  Found the following (possibly) invalid DOIs:
    DOI: 10.3102/1076998614547577
      From: DESCRIPTION
      Status: Service Unavailable
      Message: 503
      
  These URLs and DOIs are valid.

## revdepcheck results

We checked 3 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

* We saw 0 new problems
* We failed to check 0 packages
