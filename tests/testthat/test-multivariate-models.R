library(dplyr, warn.conflicts = FALSE, quietly = TRUE)
library(tidyr, quietly = TRUE)
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

# mod <- bdf_MVML
# struct <- mod$modelStruct$corStruct

test_that("targetVariance() works with multivariate models.", {
  test_Sigma_mats(bdf_MVML, bdf_long$schoolNR)
})

test_that("Derivative matrices are of correct dimension with multivariate models.", {
  test_deriv_dims(bdf_MVML)
})

test_that("Information matrices work with FIML too.", {
  test_with_FIML(bdf_MVML)
})
