## Resubmission 2
This resubmission omits the space within the doi specification in Description.

## Resubmission
This is a resubmission. This resubmission adds undirected single quotes for non-English usage and package names in title and description in the DESCRIPTION file. The resubmission also fixes the abbreviations and extra spaces between words in description in the DESCRIPTION file.

## Test environments
* local Windows 10 Pro, R 3.6.3
* Ubuntu 16.04.6 LTS (on travis-ci), R 3.6.2
* win-builder (devel, release, oldrelease)

## R CMD check results
There were no ERRORs or WARNINGs.

There were two NOTEs:
  
* Possibly mis-spelled words in DESCRIPTION: Pustejovsky (24:48) Shadish (24:73) gls (16:77, 25:28) lme (16:34, 20:57, 25:21)
  
  All of the identified words are spelled correctly.

* Found the following (possibly) invalid URLs: URL: https://doi.org/10.2307/2533274 
  (moved to https://www.jstor.org/stable/2533274) 
  From: inst/doc/Information-matrices-for-fitted-LME-models.html README.md 
  Status: 403 Message: Forbidden
  
  The flagged URL is correct.
  
