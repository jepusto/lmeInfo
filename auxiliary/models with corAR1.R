library(nlme)

data(Laski, package = "scdhlm")

Laski_AR1 <- lme(fixed = outcome ~ treatment,
                 random = ~ treatment | case,
                 correlation = corAR1(0.2, ~ time | case),
                 data = Laski)

inherits(Laski_AR1$modelStruct$corStruct, "corAR1") # TRUE

Laski_hetAR1 <- lme(fixed = outcome ~ treatment,
                    random = ~ treatment | case,
                    correlation = corAR1(0, ~ time | case),
                    weights = varIdent(form = ~ 1 | treatment),
                    data = Laski)

inherits(Laski_hetAR1$modelStruct$corStruct, "corAR1") # TRUE


data(egsingle, package = "mlmRev")

lme_2level <-
  lme(fixed = math ~ year + female + black + hispanic,
      random = ~ 1 | childid,
      correlation = corAR1(0, ~ year | childid),
      data = egsingle)

inherits(lme_2level$modelStruct$corStruct, "corAR1") # FALSE


lme_3level <-
  lme(fixed = math ~ year + female + black + hispanic,
      random = ~ 1 | schoolid/childid,
      correlation = corAR1(0.1, ~ year | schoolid/childid),
      data = egsingle)

inherits(lme_3level$modelStruct$corStruct, "corAR1") # FALSE
