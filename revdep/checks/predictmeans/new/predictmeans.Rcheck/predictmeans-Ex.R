pkgname <- "predictmeans"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('predictmeans')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CookD")
### * CookD

flush(stderr()); flush(stdout())

### Name: CookD
### Title: Calculates and plots Cook's distances for a Linear (Mixed) Model
### Aliases: CookD

### ** Examples

library(predictmeans)
Oats$nitro <- factor(Oats$nitro)
fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# library(lme4)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
CookD(fm)



cleanEx()
nameEx("Kmatrix")
### * Kmatrix

flush(stderr()); flush(stdout())

### Name: Kmatrix
### Title: Matrix of Coefficients in a Linear Model
### Aliases: Kmatrix

### ** Examples

  library(predictmeans)
  data(Oats, package="nlme")
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
  fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
  Kmatrix(fm, "Variety", prtnum=TRUE)$K
  Kmatrix(fm, "Variety", 0.5, prtnum=TRUE)$K
 # Kmatrix(fm, "Variety", "nitro")$K
  Kmatrix(fm, "Variety", "nitro", covariateV=seq(0, 0.6, 0.1))$K



cleanEx()
nameEx("PMplot")
### * PMplot

flush(stderr()); flush(stdout())

### Name: PMplot
### Title: Level Plot of a Matrix of p-values.
### Aliases: PMplot

### ** Examples

  library(predictmeans)
  set.seed(2013)
  pvalues <- runif(28)
  pmatrix <- matrix(0,8,8)
  pmatrix[lower.tri(pmatrix)] <- pvalues
  round(pmatrix, 4)
  PMplot(pmatrix)

  Oats$nitro <- factor(Oats$nitro)
  fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
  predictout <- predictmeans(fm, "nitro:Variety", atvar="Variety", adj="BH", barplot=TRUE)
  PMplot(predictout$p_valueMatrix)
  



cleanEx()
nameEx("anovalmer")
### * anovalmer

flush(stderr()); flush(stdout())

### Name: anovalmer
### Title: ANOVA of a Linear Mixed Effects Model produced by 'lmer'
###   function
### Aliases: anovalmer

### ** Examples

## Not run for simplifying process of submiting pkg to CRAN
library(predictmeans)
Oats$nitro <- factor(Oats$nitro) 
fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
anovalmer(fm)



cleanEx()
nameEx("contrastmeans")
### * contrastmeans

flush(stderr()); flush(stdout())

### Name: contrastmeans
### Title: Linear Contrast Tests for a Linear Model
### Aliases: contrastmeans

### ** Examples

library(predictmeans)
# ftable(xtabs(yield ~ Block+Variety+nitro, data=Oats))
Oats$nitro <- factor(Oats$nitro)
fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# library(lme4)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)

## Not run: 
## The contrast has a contrast matrix as follows:
#     0:Golden Rain 0:Marvellous 0:Victory 
#[1,]            -1            0         1 
#[2,]             0            0         1 
#     0.2:Golden Rain 0.2:Marvellous 0.2:Victory 
#[1,]               0              0           0 
#[2,]               0              0           0 
#     0.4:Golden Rain  0.4:Marvellous 0.4:Victory
#[1,]               0               0           0
#[2,]               0              -1           0
#      0.6:Golden Rain 0.6:Marvellous 0.6:Victory
#[1,]                0              0           0
#[2,]                0              0           0

# 1. Enter above contrast matrix into a pop up window, then close the window
# contrastmeans(fm, "nitro:Variety")
 
# 2. Construct the contrast matrix directly
cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
            c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
contrastmeans(fm, "nitro:Variety", ctrmatrix=cm)



cleanEx()
nameEx("covariatemeans")
### * covariatemeans

flush(stderr()); flush(stdout())

### Name: covariatemeans
### Title: Predicted Means of a Linear Model with Covariate Variable(s)
### Aliases: covariatemeans

### ** Examples

  library(predictmeans)
  data(Oats, package="nlme")
  fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# library(lme4)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
  covariatemeans(fm, "Variety", covariate="nitro")
  covariatemeans(fm, "Variety", covariate="nitro", covariateV=seq(0, 0.6, 0.1))$data



cleanEx()
nameEx("permanova.lmer")
### * permanova.lmer

flush(stderr()); flush(stdout())

### Name: permanova.lmer
### Title: Permutation ANOVA for 'lmer' Model
### Aliases: permanova.lmer

### ** Examples

library(predictmeans)
Oats$nitro <- factor(Oats$nitro) 
fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)

## Permutation Test for model terms
# permanova.lmer(fm)
# permanova.lmer(fm, drop=FALSE)
## Compare to F test
# fm0 <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# anova(fm0)



cleanEx()
nameEx("permindex")
### * permindex

flush(stderr()); flush(stdout())

### Name: permindex
### Title: Permutation Index
### Aliases: permindex

### ** Examples

  library(predictmeans)
  block <- rep(1:3, each=12)
  group <- rep(rep(1:3, each=4), 3)
  data <- data.frame(block, group)
  cbind(data, permindex(data, block="block", group="group", nsim=5))  
                        # Permute group as a whole within each block first, 
                        # then permute obs within each group.
  cbind(data, permindex(data, block="block",  nsim=5)) 
                        # Permute obs within each block only.
  cbind(data, permindex(data, group="group", nsim=5)) 
                        # Permute groups as a whole block first, 
                        # then permute obs within each group.
  cbind(data, permindex(data, nsim=5))  # Free permutation.



cleanEx()
nameEx("permlmer")
### * permlmer

flush(stderr()); flush(stdout())

### Name: permlmer
### Title: Permutation Test of random or fixed effects for 'lmer' model.
### Aliases: permlmer

### ** Examples

# library(predictmeans)
## Test random effects
# fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
# fm2 <- lmer(Reaction ~ Days + (Days || Subject), sleepstudy)
# fm3 <- update(fm1, . ~ . - (Days | Subject) + (1 | Subject))
# anova(fm1, fm2, fm3)
# permlmer(fm3, fm2)
# permlmer(fm2, fm1)

## Test fixed effects
# Oats$nitro <- factor(Oats$nitro)
# fm0 <- lmer(yield ~ nitro+Variety+(1|Block/Variety), data=Oats)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
# permlmer(fm0, fm)



cleanEx()
nameEx("permmodels")
### * permmodels

flush(stderr()); flush(stdout())

### Name: permmodels
### Title: Permutation Test of Linear Model
### Aliases: permmodels

### ** Examples

## Not run for simplifying process of submiting pkg to CRAN
#library(predictmeans)
#Oats$nitro <- factor(Oats$nitro) 
#fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
## library(lme4)
## fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#
## Permutation Test for model terms
#system.time(
#  permlme <- permmodels(model=fm, data=Oats, block="Block", group="Variety", nsim=999)
#)  
#
## Permutation Test for multiple comparisons
#predictmeans(model=fm, modelterm="nitro:Variety", atvar="Variety", adj="BH", 
#  permlist=permlme, plot=FALSE)
#
## Permutation Test for specified contrasts
#cm <- rbind(c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
#            c(0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0))
#contrastmeans(model=fm, modelterm="nitro:Variety", ctrmatrix=cm, permlist=permlme)



cleanEx()
nameEx("predictmeans")
### * predictmeans

flush(stderr()); flush(stdout())

### Name: predictmeans
### Title: Predicted Means of a Linear Model
### Aliases: predictmeans

### ** Examples

  library(predictmeans)
  ftable(xtabs(yield ~ Block+Variety+nitro, data=Oats))
  Oats$nitro <- factor(Oats$nitro)
  fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# library(lme4)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
  predictmeans(fm, "nitro", adj="BH")
  predictmeans(fm, "nitro:Variety", atvar="Variety", adj="BH")
  predictout <- predictmeans(fm, "nitro:Variety", atvar="Variety", adj="BH", barplot=TRUE)
  names(predictout)
  print(predictout$predictmeansPlot)
  print(predictout$predictmeansBarPlot)



cleanEx()
nameEx("residplot")
### * residplot

flush(stderr()); flush(stdout())

### Name: residplot
### Title: Diagnostic Plots for a Linear (Mixed) Model
### Aliases: residplot

### ** Examples

## Note that the order of levels of nested random effects is oposite 
## between lme and lmer objects.

library(predictmeans)
Oats$nitro <- factor(Oats$nitro)
fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
residplot(fm, level=2)    #lme: level=2 for random effect "Block:Variety"

#  Not Run
#  library(lme4)
#  fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
#  residplot(fm) # lmer: By default level=1 for random effect "Block:Variety"



cleanEx()
nameEx("varcomp")
### * varcomp

flush(stderr()); flush(stdout())

### Name: varcomp
### Title: Calculate stder and CI of variance components for 'lmer' or
###   'lme' model
### Aliases: varcomp

### ** Examples

library(predictmeans)
Oats$nitro <- factor(Oats$nitro) 
fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
## Not run: varcomp(fm)
fm1 <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
varcomp(fm1)

data(Orthodont, package="nlme")
mod <- lmer(distance ~ age + (age|Subject), data=Orthodont)
## Not run: varcomp(mod)
mod1 <- lme(distance ~ age, random=~age|Subject, data=Orthodont)
varcomp(mod1)



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
