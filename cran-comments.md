## Test environments
* local Windows 10 Pro, R 3.6.3
* Ubuntu 16.04.6 LTS (on travis-ci), R 3.6.2
* win-builder (devel, release, oldrelease)

## R CMD check results
There were no ERRORs or WARNINGs.

There were two NOTEs:
  
* Possibly mis-spelled words in DESCRIPTION: Pustejovsky (23:35) Shadish (23:60) lme (20:10) lmeStruct (3:33)
  
  All of the identified words are spelled correctly.

* Found the following (possibly) invalid URLs: URL: https://doi.org/10.2307/2533274 
  (moved to https://www.jstor.org/stable/2533274) 
  From: inst/doc/Information-matrices-for-fitted-LME-models.html README.md 
  Status: 403 Message: Forbidden
  
  The flagged URL is correct.
  
