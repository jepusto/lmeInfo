library(nlme, quietly=TRUE, warn.conflicts=FALSE)

data(Hartnagel, package = "carData")

Hart_AR <- gls(fconvict ~ tfr + partic + degrees + mconvict,
                correlation=corAR1(0.3), method="REML",
                data=Hartnagel)
Hart_CAR <- gls(fconvict ~ tfr + partic + degrees + mconvict,
                correlation=corAR1(0.3), method="REML",
                data=Hartnagel)

Hart_AR1 <- gls(fconvict ~ tfr + partic + degrees + mconvict,
                correlation=corAR1(0.3, form = ~ year), method="REML",
                data=Hartnagel)
Hart_CAR1 <- gls(fconvict ~ tfr + partic + degrees + mconvict,
                 correlation=corCAR1(0.3, form = ~ year), method="REML",
                 data=Hartnagel)
Hart_MA1 <- gls(fconvict ~ tfr + partic + degrees + mconvict,
                correlation=corARMA(0.3, form = ~ year, p=0,q=1), method="REML",
                data=Hartnagel)


data(Orthodont)

Ortho_hom <- gls(distance ~ Subject:age + Sex,
                 data = Orthodont)

Ortho_power <- gls(distance ~ Subject:age + Sex,
                   weights = varPower(),
                   data = Orthodont)

Ortho_AR1 <- gls(distance ~ age + Sex,
                 correlation = corAR1(0.2, form = ~ age | Subject),
                 data = Orthodont)

Ortho_AR1_power <- gls(distance ~ age + Sex,
                       correlation = corAR1(0.2, form = ~ age | Subject),
                       weights = varPower(),
                       data = Orthodont)

Ortho_CAR1 <- gls(distance ~ age + Sex,
                 correlation = corCAR1(0.2, form = ~ age | Subject),
                 data = Orthodont)

Ortho_CAR1_power <- gls(distance ~ age + Sex,
                       correlation = corCAR1(0.2, form = ~ age | Subject),
                       weights = varPower(),
                       data = Orthodont)

data(Laski, package = "scdhlm")

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

mod <- Hart_AR1
grps <- rep("A",nrow(Hartnagel))
invert <- TRUE
sigma_scale <- TRUE
R_list <- build_corr_mats(mod)

test_that("targetVariance() works with gls models.", {
  test_Sigma_mats(Hart_AR, rep("A", nrow(Hartnagel)))
  test_Sigma_mats(Hart_CAR, rep("A", nrow(Hartnagel)))
  test_Sigma_mats(Hart_AR1, rep("A", nrow(Hartnagel)))
  test_Sigma_mats(Hart_CAR1, rep("A", nrow(Hartnagel)))
  test_Sigma_mats(Hart_MA1, rep("A", nrow(Hartnagel)))
  test_Sigma_mats(Ortho_hom, 1:nrow(Orthodont))
  test_Sigma_mats(Ortho_power, 1:nrow(Orthodont))
  test_Sigma_mats(Ortho_AR1, Orthodont$Subject)
  test_Sigma_mats(Ortho_AR1_power, Orthodont$Subject)
  test_Sigma_mats(Ortho_CAR1, Orthodont$Subject)
  test_Sigma_mats(Ortho_CAR1_power, Orthodont$Subject)
  test_Sigma_mats(Laski_AR1, Laski$case)
  test_Sigma_mats(Laski_het, 1:nrow(Laski))
  test_Sigma_mats(Laski_hetAR1, Laski$case)
  test_Sigma_mats(Laski_CAR1, Laski$case)
})

test_that("Derivative matrices are of correct dimension with gls models.", {
  test_deriv_dims(Hart_AR)
  test_deriv_dims(Hart_CAR)
  test_deriv_dims(Hart_AR1)
  test_deriv_dims(Hart_CAR1)
  test_deriv_dims(Hart_MA1)
  test_deriv_dims(Ortho_hom)
  test_deriv_dims(Ortho_power)
  test_deriv_dims(Ortho_AR1)
  test_deriv_dims(Ortho_AR1_power)
  test_deriv_dims(Ortho_CAR1)
  test_deriv_dims(Ortho_CAR1_power)
  test_deriv_dims(Laski_AR1)
  test_deriv_dims(Laski_het)
  test_deriv_dims(Laski_hetAR1)
  test_deriv_dims(Laski_CAR1)
})

test_that("Information matrices work with FIML too.", {
  test_with_FIML(Hart_AR)
  test_with_FIML(Hart_CAR)
  test_with_FIML(Hart_AR1)
  test_with_FIML(Hart_CAR1)
  test_with_FIML(Hart_MA1)
  test_with_FIML(Ortho_hom)
  test_with_FIML(Ortho_power)
  test_with_FIML(Ortho_AR1)
  test_with_FIML(Ortho_AR1_power)
  test_with_FIML(Ortho_CAR1)
  test_with_FIML(Ortho_CAR1_power)
  test_with_FIML(Laski_AR1)
  test_with_FIML(Laski_het)
  test_with_FIML(Laski_hetAR1)
  test_with_FIML(Laski_CAR1)
})

test_that("dR_dcorStruct.corCAR1 returns the same result as dR_dcorStruct.corAR1.", {
  expect_equal(dR_dcorStruct.corCAR1(Hart_AR1$modelStruct$corStruct),
               dR_dcorStruct.corAR1(Hart_CAR1$modelStruct$corStruct),
               tol = 10^-5)
  expect_equal(dR_dcorStruct.corCAR1(Hart_AR$modelStruct$corStruct),
               dR_dcorStruct.corAR1(Hart_CAR$modelStruct$corStruct),
               tol = 10^-5)
  expect_equal(dR_dcorStruct.corCAR1(Ortho_AR1$modelStruct$corStruct),
               dR_dcorStruct.corAR1(Ortho_CAR1$modelStruct$corStruct),
               tol = 10^-5)
  expect_equal(dR_dcorStruct.corCAR1(Ortho_AR1_power$modelStruct$corStruct),
               dR_dcorStruct.corAR1(Ortho_CAR1_power$modelStruct$corStruct),
               tol = 10^-5)
  expect_equal(dR_dcorStruct.corCAR1(Laski_AR1$modelStruct$corStruct),
               dR_dcorStruct.corAR1(Laski_CAR1$modelStruct$corStruct),
               tol = 10^-4)

})


test_that("Results do not depend on order of data.", {
  test_after_shuffling(Hart_AR, seed = 20)
  test_after_shuffling(Hart_CAR, seed = 20)
  test_after_shuffling(Hart_AR1, seed = 20)
  test_after_shuffling(Hart_CAR1, seed = 20)
  test_after_shuffling(Hart_MA1, seed = 20)
  test_after_shuffling(Ortho_hom, seed = 20)
  test_after_shuffling(Ortho_power, seed = 20)
  test_after_shuffling(Ortho_AR1, seed = 20)
  test_after_shuffling(Ortho_AR1_power, seed = 20)
  test_after_shuffling(Ortho_CAR1, seed = 20)
  test_after_shuffling(Ortho_CAR1_power, seed = 20)
  test_after_shuffling(Laski_AR1, seed = 20)
  test_after_shuffling(Laski_het, seed = 21)
  test_after_shuffling(Laski_hetAR1, seed = 20)
  test_after_shuffling(Laski_CAR1, seed = 20)
})


test_that("New REML calculations work.", {

  check_REML2(Hart_AR)
  check_REML2(Hart_CAR)
  check_REML2(Hart_AR1)
  check_REML2(Hart_CAR1)
  check_REML2(Hart_MA1)
  check_REML2(Ortho_hom)
  check_REML2(Ortho_power)
  check_REML2(Ortho_AR1)
  check_REML2(Ortho_AR1_power)
  check_REML2(Ortho_CAR1)
  check_REML2(Ortho_CAR1_power)
  check_REML2(Laski_AR1)
  check_REML2(Laski_het)
  check_REML2(Laski_hetAR1)
  check_REML2(Laski_CAR1)
})
