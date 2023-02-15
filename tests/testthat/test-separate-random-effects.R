library(nlme)

skip_if_not_installed("scdhlm")

data(Bryant2018, package = "scdhlm")
Bryant2018 <- subset(Bryant2018, !is.na(outcome))

#------------------------------------------------------------------------------
# varying intercept at level 3, ALLOW intercept & trend covariance at level 2
#------------------------------------------------------------------------------

# One way to specify the random effects
mod1_1 <- lme(fixed = outcome ~ session_c + treatment + session_trt,
              random = list( ~ 1 | school, ~ session_c | case),
              correlation = corAR1(0, ~ session_c | school/case),
              data = Bryant2018)
g_mod1_1 <- g_mlm(mod1_1, p_const = c(0,0,1,17), r_const = c(1,1,0,0,0,1))

# Second way: names list
# the names define the grouping factors and the formulas describe the random-effects models at each level
mod1_2 <- lme(fixed = outcome ~ session_c + treatment + session_trt,
              random = list(school = ~ 1, case = ~ session_c),
              correlation = corAR1(0, ~ session_c | school/case),
              data = Bryant2018)
g_mod1_2 <- g_mlm(mod1_2, p_const = c(0,0,1,17), r_const = c(1,1,0,0,0,1))

test_that("mod1_1 and mod1_2 returns the same results.", {
  expect_equal(mod1_1$coefficients$fixed, mod1_2$coefficients$fixed)
  expect_equal(mod1_1$coefficients$random, mod1_2$coefficients$random)
  expect_equal(mod1_1$sigma, mod1_2$sigma)
  expect_equal(mod1_1$varFix, mod1_2$varFix)
  expect_equal(mod1_1$logLik, mod1_2$logLik)
  expect_equal(g_mod1_1$g_AB, g_mod1_2$g_AB)
  expect_equal(g_mod1_1$SE_g_AB, g_mod1_2$SE_g_AB)
})

#-------------------------------------------------------------------------------------
# varying intercept at level 3, DO NOT ALLOW intercept & trend covariance at level 2
#-------------------------------------------------------------------------------------

# first approach to specify the random effects
mod2_1 <- suppressWarnings(lme(fixed = outcome ~ session_c + treatment + session_trt,
              random = list( ~ 1 | school, ~ 1 | case, ~ 0 + session_c | case),
              correlation = corAR1(0, ~ session_c | school/case),
              data = Bryant2018))
# warning: cannot use smaller level of grouping for 'correlation' than for 'random'. Replacing the former with the latter.
VarCorr(mod2_1) # pdLogChol(1) parametrization

# Fisher_info(mod2_1)
# g_mod2_1 <- g_mlm(mod2_1, p_const = c(0,0,1,17), r_const = c(1,0,1,0,1))


# second approach: pdDiag in pdMat classes
mod2_2 <- lme(fixed = outcome ~ session_c + treatment + session_trt,
              random = list(school = ~ 1, case = pdDiag(~ session_c)),
              correlation = corAR1(0, ~ session_c | school/case),
              data = Bryant2018)
VarCorr(mod2_2) # pdDiag(1) parametrization?

Fisher_info(mod2_2)
g_mod2_2 <- g_mlm(mod2_2, p_const = c(0,0,1,17), r_const = c(1,1,0,0,1))
#summary(g_mod2_2)

test_that("mod2_1 and mod2_2 return the same results.", {
  expect_equal(mod2_1$coefficients$fixed, mod2_2$coefficients$fixed, tol = 1e-5)
  expect_equal(as.numeric(VarCorr(mod2_1)[,1][2]), as.numeric(VarCorr(mod2_2)[,1][2]), tol = 1e-3)
  expect_equal(as.numeric(VarCorr(mod2_1)[,2][4]), as.numeric(VarCorr(mod2_2)[,2][4]), tol = 1e-7)
  expect_equal(as.numeric(VarCorr(mod2_1)[,2][6]), as.numeric(VarCorr(mod2_2)[,2][5]), tol = 1e-7)
  expect_equal(as.numeric(VarCorr(mod2_1)[,2][7]), as.numeric(VarCorr(mod2_2)[,2][6]))
  expect_equal(mod2_1$sigma, mod2_2$sigma)
  expect_equal(mod2_1$varFix, mod2_2$varFix, tolerance = 5e-4) # intercept var differs at 6th decimals
  expect_equal(mod2_1$logLik, mod2_2$logLik)
})


test_that("Fisher_info() of mod1_1 and mod2_2 return similar output.", {
  info1 <- Fisher_info(mod1_1)[-3, -3] # remove the row and col with the cor btw intercept and slope involved
  info2 <- Fisher_info(mod2_2)
  expect_equivalent(info1 > 0, info2 > 0)
  expect_equal(sum(diag(info1)), sum(diag(info2)), tol = 100)
  expect_equal(det(info1), det(info2), tol = 1e-3)
  # expect_equal(info1[, c(1:2, 5)], info2[, c(1:2, 5)], tolerance = 4e-2)
})


test_that("targetVariance() works with separate random effects models.", {
  test_Sigma_mats(mod1_1, Bryant2018$school)
  test_Sigma_mats(mod1_2, Bryant2018$school)
  test_Sigma_mats(mod2_2, Bryant2018$school)
})


test_that("Derivative matrices are of correct dimension with separate random effects models.", {
  test_deriv_dims(mod1_1)
  test_deriv_dims(mod1_2)
  test_deriv_dims(mod2_2)
})


test_that("Information matrices work with FIML with separate random effects models.", {
  test_with_FIML(mod1_1)
  test_with_FIML(mod1_2)
  test_with_FIML(mod2_2)
})

test_that("Info matrices work with dropped observations.", {

  skip_on_cran()

  test_after_deleting(mod1_1)
  test_after_deleting(mod1_2)
  test_after_deleting(mod2_2)

})
