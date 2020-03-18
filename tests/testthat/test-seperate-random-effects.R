library(nlme)

data(Bryant2018)

#------------------------------------------------------------------------------
# varying intercept at level 3, ALLOW intercept & trend covariance at level 2
#------------------------------------------------------------------------------

# One way to specify the random effects
mod1_1 <- lme(fixed = outcome ~ session_c + treatment + session_trt,
              random = list( ~ 1 | school, ~ session_c | case),
              correlation = corAR1(0, ~ session_c | school/case),
              data = Bryant2018)
summary(mod1_1)
g_mod1_1 <- g_mlm(mod1_1, p_const = c(0,0,1,17), r_const = c(1,1,0,0,0,1))
#print(g_mod1_1)

# Second way: names list
# the names define the grouping factors and the formulas describe the random-effects models at each level
mod1_2 <- lme(fixed = outcome ~ session_c + treatment + session_trt,
              random = list(school = ~ 1, case = ~ session_c),
              correlation = corAR1(0, ~ session_c | school/case),
              data = Bryant2018)
summary(mod1_2)
g_mod1_2 <- g_mlm(mod1_2, p_const = c(0,0,1,17), r_const = c(1,1,0,0,0,1))
#print(g_mod1_2)

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
mod2_1 <- lme(fixed = outcome ~ session_c + treatment + session_trt,
              random = list( ~ 1 | school, ~ 1 | case, ~ 0 + session_c | case),
              correlation = corAR1(0, ~ session_c | school/case),
              data = Bryant2018)
# warning: cannot use smaller level of grouping for 'correlation' than for 'random'. Replacing the former with the latter.
summary(mod2_1)
VarCorr(mod2_1) # pdLogChol(1) parametrization
# intervals(mod2_1)

Fisher_info(mod2_1)
g_mod2_1 <- g_mlm(mod2_1, p_const = c(0,0,1,17), r_const = c(1,0,1,0,1))
# issue: three groups (1 for school level, 2 for case level)
summary(g_mod2_1)
print(g_mod2_1)


# second approach: pdDiag in pdMat classes
mod2_2 <- lme(fixed = outcome ~ session_c + treatment + session_trt,
              random = list(school = ~ 1, case = pdDiag(~ session_c)),
              correlation = corAR1(0, ~ session_c | school/case),
              data = Bryant2018)
summary(mod2_2)
VarCorr(mod2_2) # pdDiag(1) parametrization?

Fisher_info(mod2_2)
g_mod2_2 <- g_mlm(mod2_2, p_const = c(0,0,1,17), r_const = c(1,1,0,0,1))
summary(g_mod2_2)
print(g_mod2_2)

test_that("mod1_1 and mod1_2 returns the same results.", {
  expect_equal(mod2_1$coefficients$fixed, mod2_2$coefficients$fixed)
  expect_equal(as.numeric(VarCorr(mod2_1)[,2][2]), as.numeric(VarCorr(mod2_2)[,2][2]), tolerance = 1e-4)
  expect_equal(as.numeric(VarCorr(mod2_1)[,2][4]), as.numeric(VarCorr(mod2_2)[,2][4]), tolerance = 1e-7)
  expect_equal(as.numeric(VarCorr(mod2_1)[,2][6]), as.numeric(VarCorr(mod2_2)[,2][5]), tolerance = 1e-7)
  expect_equal(as.numeric(VarCorr(mod2_1)[,2][7]), as.numeric(VarCorr(mod2_2)[,2][6]))
  expect_equal(mod2_1$sigma, mod2_2$sigma)
  expect_equal(mod2_1$varFix, mod2_2$varFix, tolerance = 1e-7) # intercept var differs at 6th decimals
  expect_equal(mod2_1$logLik, mod2_2$logLik)
})


# third approach
mod2_3 <- lme(fixed = outcome ~ session_c + treatment + session_trt,
              random = list( ~ 1 | school, ~ 1 | case, ~ 0 + session_c | case),
              correlation = corAR1(0, ~ session_c | school/case/case), # get rid of the warning in mod2_1
              data = Bryant2018)
summary(mod2_3)
