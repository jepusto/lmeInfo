pkgname <- "merDeriv"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('merDeriv')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("bread.glmerMod")
### * bread.glmerMod

flush(stderr()); flush(stdout())

### Name: bread.glmerMod
### Title: Extract Bread Component for Huber-White Sandwich Estimator of
###   Generalized Linear Mixed Effects Models
### Aliases: bread.glmerMod

### ** Examples

## Not run: 
##D # The cbpp example
##D data(finance, package = "smdata")
##D 
##D lme4fit <- glmer(corr ~ jmeth + (1 | item), data = finance,
##D                  family = binomial, nAGQ = 20)
##D 
##D # bread component for all parameters
##D bread(lme4fit, full = TRUE, ranpar = "var")
## End(Not run)



cleanEx()
nameEx("bread.lmerMod")
### * bread.lmerMod

flush(stderr()); flush(stdout())

### Name: bread.lmerMod
### Title: Extract Bread Component for Huber-White Sandwich Estimator of
###   Linear Mixed Effects Models
### Aliases: bread.lmerMod

### ** Examples

## Not run: 
##D # The sleepstudy example
##D lme4fit <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE)
##D 
##D # bread component for all parameters
##D bread(lme4fit, full = TRUE, information = "expected", ranpar = "var")
## End(Not run)



cleanEx()
nameEx("estfun.glmerMod")
### * estfun.glmerMod

flush(stderr()); flush(stdout())

### Name: estfun.glmerMod
### Title: Extract Cluster-wise Derivatives for Generalized Linear Mixed
###   Effects Models
### Aliases: estfun.glmerMod

### ** Examples

## Not run: 
##D data(finance, package = "smdata")
##D 
##D lme4fit <- glmer(corr ~ jmeth + (1 | item), data = finance,
##D                  family = binomial, nAGQ = 20)
##D 
##D # clusterwise scores
##D estfun(lme4fit, ranpar = "var")
## End(Not run)



cleanEx()
nameEx("estfun.lmerMod")
### * estfun.lmerMod

flush(stderr()); flush(stdout())

### Name: estfun.lmerMod
### Title: Extract Case-wise and Cluster-wise Derivatives for Linear Mixed
###   Effects Models
### Aliases: estfun.lmerMod

### ** Examples

## Not run: 
##D # The sleepstudy example
##D lme4fit <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE)
##D 
##D # casewise scores
##D estfun(lme4fit, level = 1, ranpar = "var")
##D 
##D # clusterwise scores
##D estfun(lme4fit, level = 2, ranpar = "sd")
## End(Not run)



cleanEx()
nameEx("llcont.glmerMod")
### * llcont.glmerMod

flush(stderr()); flush(stdout())

### Name: llcont.glmerMod
### Title: Extract Cluster-wise Log Likelihoods for Generalized Linear
###   Mixed Effects Models
### Aliases: llcont.glmerMod

### ** Examples

## Not run: 
##D data(finance, package="smdata")
##D 
##D lme4fit <- glmer(corr ~ jmeth + (1 | item), data = finance,
##D                  family = binomial, nAGQ = 20)
##D 
##D # clusterwise log likelihood
##D llcont(lme4fit)
## End(Not run)  



cleanEx()
nameEx("llcont.lmerMod")
### * llcont.lmerMod

flush(stderr()); flush(stdout())

### Name: llcont.lmerMod
### Title: Extract Case-wise Log Likelihoods for Linear Mixed Effects
###   Models
### Aliases: llcont.lmerMod

### ** Examples

## Not run: 
##D # The sleepstudy example
##D lme4fit <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE)
##D 
##D # clusterwise log likelihood
##D llcont(lme4fit)
## End(Not run)  



cleanEx()
nameEx("vcov.glmerMod")
### * vcov.glmerMod

flush(stderr()); flush(stdout())

### Name: vcov.glmerMod
### Title: Extract Variance-Covariance Matrix of all Parameters for
###   Generalized Linear Mixed Effects Models
### Aliases: vcov.glmerMod

### ** Examples

## Not run: 
##D # The cbpp example
##D data(finance, package="smdata")
##D 
##D lme4fit <- glmer(corr ~ jmeth + (1 | item), data = finance,
##D                  family = binomial, nAGQ = 20)
##D 
##D # variance covariance matrix for all parameters
##D vcov(lme4fit, full = TRUE, ranpar = "var")
## End(Not run)



cleanEx()
nameEx("vcov.lmerMod")
### * vcov.lmerMod

flush(stderr()); flush(stdout())

### Name: vcov.lmerMod
### Title: Extract Variance-Covariance Matrix of all Parameters for Linear
###   Mixed Effects Models
### Aliases: vcov.lmerMod

### ** Examples

## Not run: 
##D # The sleepstudy example
##D lme4fit <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE)
##D 
##D # variance covariance matrix for all parameters
##D vcov(lme4fit, full = TRUE, ranpar = "var")
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
