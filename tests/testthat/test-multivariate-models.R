library(dplyr, warn.conflicts = FALSE, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(nlme)

data(bdf, package = "mlmRev")

bdf_long <-
  bdf %>%
  filter(schoolNR %in% levels(schoolNR)[1:12]) %>%
  droplevels() %>%
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


# introduce random missing to bdf_long

bdf_long_wm <-
  bdf_long %>%
  mutate(
    row_index = rbinom(n = n(), size = 1, prob = 0.9),
    score = if_else(row_index == 1, score, as.numeric(NA)),
    measure_id = as.integer(factor(measure))
  ) %>%
  filter(!is.na(score)) %>%
  select(-row_index)

bdf_wm <- lme(score ~ 0 + measure,
              random = ~ 1| schoolNR / pupilNR,
              corr = corSymm(form = ~ measure_id | schoolNR / pupilNR),
              weights = varIdent(form = ~ 1 | measure),
              data = bdf_long_wm,
              control=lmeControl(msMaxIter = 100, apVar = FALSE, returnObject = TRUE))

mod <- bdf_wm
struct <- mod$modelStruct$corStruct

test_that("targetVariance() works with multivariate models.", {
  test_Sigma_mats(bdf_wm, bdf_long_wm$schoolNR)
})

test_that("Derivative matrices are of correct dimension with multivariate models.", {
  test_deriv_dims(bdf_wm)
})

test_that("Information matrices work with FIML too.", {
  test_with_FIML(bdf_wm)
})
