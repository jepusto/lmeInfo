library(dplyr)
library(tidyr)
library(nlme)
data(bdf, package = "mlmRev")

bdf_long <-
  bdf %>%
  pivot_longer(cols = c(IQ.verb, IQ.perf, aritPRET),
               names_to = "measure",
               values_to = "score") %>%
  select(schoolNR, pupilNR, sex, Minority, measure, score)

bdf_MVML <- lme(score ~ 0 + measure,
                random = ~ 1| schoolNR / pupilNR,
                corr = corSymm(form = ~ 1 | schoolNR / pupilNR),
                weights = varIdent(form = ~ 1 | measure),
                data = bdf_long)

summary(bdf_MVML)

bdf_MVML$modelStruct$corStruct
bdf_MVML$modelStruct$varStruct

mod <- bdf_MVML

struct <- mod$modelStruct$corStruct

res <- dR_dcorStruct(struct)



