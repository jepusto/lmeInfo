library(nlme)
data(Laski)
data(Thiemann2001)

# Two-level model
Laski_fx_sigma <- lme(fixed = outcome ~ treatment,
                      random = ~ treatment | case,
                      control = lmeControl(sigma = 2),
                      data = Laski)

attr(Laski_fx_sigma$modelStruct, "fixedSigma") # Indicating whether sigma is fixed or not
Fisher_info(Laski_fx_sigma, type = "averaged")

# Three-level model
Thiemann2001_fx_sigma <- lme(fixed = outcome ~ treatment,
                             random = ~ 1 | case/series,
                             correlation = corAR1(0, ~ time | case/series),
                             control = lmeControl(sigma = 1),
                             data = Thiemann2001)

attr(Thiemann2001_fx_sigma$modelStruct, "fixedSigma")
Fisher_info(Thiemann2001_fx_sigma, type = "expected")

test_that("targetVariance() works with models when sigma is fixed.", {
  test_Sigma_mats(Laski_fx_sigma, Laski$case)
  test_Sigma_mats(Thiemann2001_fx_sigma, Thiemann2001$case)
})

test_that("Derivative matrices are of correct dimension with models when sigma is fixed.", {
  expect_error(test_deriv_dims(Laski_fx_sigma))
  expect_error(test_deriv_dims(Thiemann2001_fx_sigma))
})

test_that("Information matrices work with FIML too when sigma is fixed.", {
  expect_error(test_with_FIML(Laski_fx_sigma))
  expect_error(test_with_FIML(Thiemann2001_fx_sigma))
})
