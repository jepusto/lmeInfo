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
                 correlation = corAR1(0.2, ~ time | case),
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

Laski_MA2 <- lme(fixed = outcome ~ treatment,
                     random = ~ treatment | case,
                     correlation = corARMA(c(0,0), ~ time | case, p = 0, q = 2),
                     data = Laski)

Laski_AR1MA1 <- lme(fixed = outcome ~ treatment,
                    random = ~ treatment | case,
                    correlation = corARMA(c(0,0), ~ time | case, p = 1, q = 1),
                    data = Laski)

test_that("targetVariance() works with 2-level models.", {
  test_Sigma_mats(Laski_iid, Laski$case)
  test_Sigma_mats(Laski_het, Laski$case)
  test_Sigma_mats(Laski_AR1, Laski$case)
  test_Sigma_mats(Laski_hetAR1, Laski$case)
  test_Sigma_mats(Laski_CAR1, Laski$case)
  test_Sigma_mats(Laski_MA1, Laski$case)
  test_Sigma_mats(Laski_MA2, Laski$case)
  test_Sigma_mats(Laski_AR1MA1, Laski$case)
})

test_that("Derivative matrices are of correct dimension with 2-level models.", {
  test_deriv_dims(Laski_iid)
  test_deriv_dims(Laski_het)
  test_deriv_dims(Laski_AR1)
  test_deriv_dims(Laski_hetAR1)
  test_deriv_dims(Laski_CAR1)
  test_deriv_dims(Laski_MA1)
  expect_error(test_deriv_dims(Laski_MA2))
  expect_error(test_deriv_dims(Laski_AR1MA1))
})

test_that("Information matrices work with FIML too.", {
  test_with_FIML(Laski_iid)
  test_with_FIML(Laski_het)
  test_with_FIML(Laski_AR1)
  test_with_FIML(Laski_hetAR1)
  test_with_FIML(Laski_CAR1)
  test_with_FIML(Laski_MA1)
})

test_that("dR_dcorStruct.corCAR1 returns the same result as dR_dcorStruct.corAR1.", {
  expect_equal(dR_dcorStruct.corCAR1(Laski_CAR1$modelStruct$corStruct),
               dR_dcorStruct.corAR1(Laski_AR1$modelStruct$corStruct),
               tol = 10^-5)
})

