library(nlme)

skip_if_not_installed("scdhlm")
library(scdhlm)

data(AlberMorgan, package = "scdhlm")

# Clean data

dat <- preprocess_SCD(case = case,
                      phase = condition,
                      session = session,
                      outcome = outcome,
                      design = "MB",
                      center = 5,
                      data = AlberMorgan)
dat$session.number <- dat$session
dat$trt.dummy <- dat$trt

# Fit the model
mA1 <- lme(fixed = outcome ~ 1 + session + trt + session_trt,
           random = ~ 1 + session | case,
           correlation = corAR1(0.01, ~ session | case),
           data = dat)
mA2 <- lme(fixed = outcome ~ 1 + session.number + trt + session_trt,
           random = ~ 1 + session.number | case,
           correlation = corAR1(0.01, ~ session | case),
           data = dat)
mA3 <- lme(fixed = outcome ~ 1 + I(session^1) + trt + session_trt,
           random = ~ 1 + I(session^1) | case,
           correlation = corAR1(0.01, ~ session | case),
           data = dat)
mA4 <- lme(fixed = outcome ~ 1 + I(session.number^1) + trt + session_trt,
           random = ~ 1 + I(session.number^1) | case,
           correlation = corAR1(0.01, ~ session | case),
           data = dat)
# Calculate effect size with g_mlm()
p_const <- c(0,0,1,26)
Ar_const <- c(1,62,961,0,1)

mB1 <- lme(fixed = outcome ~ 1 + session.number + trt.dummy + session_trt,
           random = ~ 1 + trt | case,
           data = dat)
mB2 <- lme(fixed = outcome ~ 1 + session.number + trt.dummy + session_trt,
           random = ~ 1 + trt.dummy | case,
           data = dat)
mB3 <- lme(fixed = outcome ~ 1 + session.number + trt.dummy + session_trt,
           random = ~ 1 + I(trt.dummy^3) | case,
           data = dat)
Br_const <- c(1,0,0,1)


test_that("extract_varcomp() works with variable names containing '.' or '^'.", {

  V_A1 <- extract_varcomp(mA1)
  V_A2 <- extract_varcomp(mA2)
  V_A3 <- extract_varcomp(mA3)
  V_A4 <- extract_varcomp(mA4)

  expect_equivalent(V_A1, V_A2)
  expect_equivalent(V_A1, V_A3)
  expect_equivalent(V_A1, V_A4)

  V_B1 <- extract_varcomp(mB1)
  V_B2 <- extract_varcomp(mB2)
  V_B3 <- extract_varcomp(mB3)

  expect_equivalent(V_B1, V_B2)
  expect_equivalent(V_B1, V_B3)

})

test_that("Fisher_info() works with variable names containing '.' or '^'.", {

  I_A1 <- Fisher_info(mA1)
  I_A2 <- Fisher_info(mA2)
  I_A3 <- Fisher_info(mA3)
  I_A4 <- Fisher_info(mA4)

  expect_equivalent(I_A1, I_A2)
  expect_equivalent(I_A1, I_A3)
  expect_equivalent(I_A1, I_A4)

  I_B1 <- Fisher_info(mB1)
  I_B2 <- Fisher_info(mB2)
  I_B3 <- Fisher_info(mB3)

  expect_equivalent(I_B1, I_B2)
  expect_equivalent(I_B1, I_B3)

})

test_that("g_mlm() works with variable names containing '.' or '^'.", {

  g_A1 <- g_mlm(mA1, p_const = p_const, r_const = Ar_const)
  g_A2 <- g_mlm(mA2, p_const = p_const, r_const = Ar_const)
  g_A3 <- g_mlm(mA3, p_const = p_const, r_const = Ar_const)
  g_A4 <- g_mlm(mA4, p_const = p_const, r_const = Ar_const)

  expect_equivalent(g_A1, g_A2)
  expect_equivalent(g_A1, g_A3)
  expect_equivalent(g_A1, g_A4)

  g_B1 <- g_mlm(mB1, p_const = p_const, r_const = Br_const)
  g_B2 <- g_mlm(mB2, p_const = p_const, r_const = Br_const)
  g_B3 <- g_mlm(mB3, p_const = p_const, r_const = Br_const)

  expect_equivalent(g_B1, g_B2)
  expect_equivalent(g_B1, g_B3)

})
