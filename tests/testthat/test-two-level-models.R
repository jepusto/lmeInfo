library(nlme)
data(Laski)

Laski_iid <- lme(fixed = outcome ~ treatment,
                 random = ~ treatment | case,
                 data = Laski)

Laski_het <- lme(fixed = outcome ~ treatment,
                 random = ~ treatment | case,
                 weights = varIdent(form = ~ 1 | treatment),
                 data = Laski)

Laski_AR1 <- lme(fixed = outcome ~ treatment,
                 random = ~ treatment | case,
                 correlation = corAR1(0, ~ time | case),
                 data = Laski)

Laski_hetAR1 <- lme(fixed = outcome ~ treatment,
                    random = ~ treatment | case,
                    correlation = corAR1(0, ~ time | case),
                    weights = varIdent(form = ~ 1 | treatment),
                    data = Laski)

Laski_CAR1 <- lme(fixed = outcome ~ treatment,
                  random = ~ treatment | case,
                  correlation = corCAR1(0.2, ~ time | case),
                  data = Laski)

Laski_MA1 <- lme(fixed = outcome ~ treatment,
                 random = ~ treatment | case,
                 correlation = corARMA(0, ~ time | case, p = 0, q = 1),
                 data = Laski)


extract_varcomp(Laski_iid)

mod <- Laski_iid
getGroups(mod)

targetVariance(mod)

test_that("targetVariance() works with 2-level models.", {
  expect_output(targetVariance(Laski_iid))
  expect_output(targetVariance(Laski_het))
  expect_output(targetVariance(Laski_AR1))
  expect_output(targetVariance(Laski_hetAR1))
  expect_output(targetVariance(Laski_CAR1))
  expect_output(targetVariance(Laski_MA1))
})

test_that("targetVariance() works with 2-level models.", {
  expect_is(targetVariance(Laski_iid), "list")
  expect_is(targetVariance(Laski_het), "list")
  expect_is(targetVariance(Laski_AR1), "list")
  expect_is(targetVariance(Laski_hetAR1), "list")
  expect_is(targetVariance(Laski_CAR1), "list")
  expect_is(targetVariance(Laski_MA1), "list")
})
