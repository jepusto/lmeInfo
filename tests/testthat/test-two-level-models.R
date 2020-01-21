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


test_that("targetVariance() works with 2-level models.", {
  test_Sigma_mats(Laski_iid, Laski$case)
  test_Sigma_mats(Laski_het, Laski$case)
  test_Sigma_mats(Laski_AR1, Laski$case)
  test_Sigma_mats(Laski_hetAR1, Laski$case)
  test_Sigma_mats(Laski_CAR1, Laski$case)
  test_Sigma_mats(Laski_MA1, Laski$case)
})

test_that("dR_dcorStruct.corCAR1 returns the same result as dR_dcorStruct.corAR1.", {
  expect_equal(dR_dcorStruct.corCAR1(Laski_CAR1$modelStruct$corStruct),
               dR_dcorStruct.corAR1(Laski_AR1$modelStruct$corStruct))
})

