library(nlme, quietly=TRUE, warn.conflicts=FALSE)

data(Ovary, package = "nlme")

Ovary$time_int <- 1:nrow(Ovary)

Ovary_hom <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), data = Ovary)

Ovary_power <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time),
                weights = varPower(),
                data = Ovary)

Ovary_AR1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time),
              correlation = corAR1(form = ~ time_int | Mare),
              data = Ovary)

Ovary_AR1_power <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time),
                    correlation = corAR1(form = ~ time_int | Mare),
                    weights = varPower(),
                    data = Ovary)

Ovary_CAR1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time),
               correlation = corCAR1(form = ~ time_int | Mare),
               data = Ovary)

Ovary_CAR1_power <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time),
                     correlation = corCAR1(form = ~ time_int | Mare),
                     weights = varPower(),
                     data = Ovary)

Laski_AR1 <- gls(outcome ~ 0 + case + case:treatment,
                 correlation = corAR1(0.2, ~ time | case),
                 data = Laski)

Laski_het <- gls(outcome ~ 0 + case + case:treatment,
                  weights = varIdent(form = ~ 1 | treatment),
                  data = Laski)

Laski_hetAR1 <- gls(outcome ~ 0 + case + case:treatment,
                    correlation = corAR1(0.2, ~ time | case),
                    weights = varIdent(form = ~ 1 | treatment),
                    data = Laski)

Laski_CAR1 <- gls(outcome ~ 0 + case + case:treatment,
                  correlation = corCAR1(0.2, ~ time | case),
                  data = Laski)

mod <- Ovary_hom
invert <- TRUE
sigma_scale <- TRUE
grps <- 1:nrow(Ovary)

test_that("targetVariance() works with gls models.", {
  test_Sigma_mats(Ovary_hom, 1:nrow(Ovary))
  test_Sigma_mats(Ovary_power, Ovary$Mare)
  test_Sigma_mats(Ovary_AR1, Ovary$Mare)
  test_Sigma_mats(Ovary_AR1_power, Ovary$Mare)
  test_Sigma_mats(Ovary_CAR1, Ovary$Mare)
  test_Sigma_mats(Ovary_CAR1_power, Ovary$Mare)
  test_Sigma_mats(Laski_AR1, Laski$case)
  test_Sigma_mats(Laski_het, Laski$case)
  test_Sigma_mats(Laski_hetAR1, Laski$case)
  test_Sigma_mats(Laski_CAR1, Laski$case)
})

test_that("Derivative matrices are of correct dimension with gls models.", {
  test_deriv_dims(Ovary_hom)
  test_deriv_dims(Ovary_power)
  test_deriv_dims(Ovary_AR1)
  test_deriv_dims(Ovary_AR1_power)
  test_deriv_dims(Ovary_CAR1)
  test_deriv_dims(Ovary_CAR1_power)
  test_deriv_dims(Laski_AR1)
  test_deriv_dims(Laski_het)
  test_deriv_dims(Laski_hetAR1)
  test_deriv_dims(Laski_CAR1)
})

test_that("Information matrices work with FIML too.", {
  test_with_FIML(Ovary_hom)
  test_with_FIML(Ovary_power)
  test_with_FIML(Ovary_AR1)
  test_with_FIML(Ovary_AR1_power)
  test_with_FIML(Ovary_CAR1)
  test_with_FIML(Ovary_CAR1_power)
  test_with_FIML(Laski_AR1)
  test_with_FIML(Laski_het)
  test_with_FIML(Laski_hetAR1)
  test_with_FIML(Laski_CAR1)
})

test_that("dR_dcorStruct.corCAR1 returns the same result as dR_dcorStruct.corAR1.", {
  expect_equal(dR_dcorStruct.corCAR1(Laski_CAR1$modelStruct$corStruct),
               dR_dcorStruct.corAR1(Laski_AR1$modelStruct$corStruct),
               tol = 10^-5)
})


test_that("Results do not depend on order of data.", {
  test_after_shuffling(Ovary_hom, seed = 20)
  test_after_shuffling(Ovary_power, seed = 20)
  test_after_shuffling(Ovary_AR1, seed = 20)
  test_after_shuffling(Ovary_AR1_power, seed = 20)
  test_after_shuffling(Ovary_CAR1, seed = 20)
  test_after_shuffling(Ovary_CAR1_power, seed = 20)
  test_after_shuffling(Laski_AR1, seed = 20)
  test_after_shuffling(Laski_het, seed = 20)
  test_after_shuffling(Laski_hetAR1, seed = 20)
  test_after_shuffling(Laski_CAR1, seed = 20)
})


test_that("New REML calculations work.", {

  check_REML2(Laski_iid)
  check_REML2(Ovary_hom)
  check_REML2(Ovary_power)
  check_REML2(Ovary_AR1)
  check_REML2(Ovary_AR1_power)
  check_REML2(Ovary_CAR1)
  check_REML2(Ovary_CAR1_power)
  check_REML2(Laski_AR1)
  check_REML2(Laski_het)
  check_REML2(Laski_hetAR1)
  check_REML2(Laski_CAR1)
})
