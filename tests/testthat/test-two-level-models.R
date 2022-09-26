library(nlme)

skip_if_not_installed("scdhlm")

data(Laski, package = "scdhlm")

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
                    correlation = corAR1(0.2, ~ time | case),
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

Laski_AR1MA1 <- lme(fixed = outcome ~ treatment,
                    random = ~ treatment | case,
                    correlation = corARMA(c(0,0), ~ time | case, p = 1, q = 1),
                    data = Laski)

Bodyweight_Gaus <- lme(weight ~ Time * Diet,
                       data = BodyWeight,
                       random = ~ Time,
                       weights = varPower(),
                       correlation = corGaus(form = ~ Time))


data(Lambert, package = "scdhlm")
Lambert_RML <- lme(fixed = outcome ~ treatment,
                   random = ~ 1 | case,
                   correlation = corAR1(0.1, ~ time | case),
                   data = subset(Lambert, measure=="academic response"))

data(Anglesea, package = "scdhlm")
Anglesea_RML <- lme(fixed = outcome ~ condition,
                    random = ~ condition | case,
                    correlation = corAR1(0.2, ~ session | case),
                    data = Anglesea)

data(Saddler, package = "scdhlm")
Saddler_quality_RML <- lme(fixed = outcome ~ treatment,
                           random = ~ 1 | case,
                           correlation = corAR1(0, ~ time | case),
                           data = subset(Saddler, measure=="writing quality"))


data(Schutte, package = "scdhlm")
Schutte <- subset(Schutte, case != "Case 4")
Schutte$case <- factor(Schutte$case)
Schutte$trt_week <- with(Schutte,
                         unlist(tapply((treatment=="treatment") * week,
                                       list(treatment, case),
                                       function(x) x - min(x))) + (treatment=="treatment"))
Schutte$week <- Schutte$week - 9 # center at follow-up time (week 9)

# Varying intercepts, fixed treatment effect, fixed trends
Schutte_RML3 <- lme(fixed = fatigue ~ week + treatment + trt_week,
                    random = ~ 1 | case,
                    correlation = corAR1(0, ~ week | case),
                    data = Schutte,
                    method = "REML")

# Varying intercepts, fixed treatment effect, varying trends
Schutte_RML4 <- update(Schutte_RML3, random = ~ week | case)

# Varying intercepts, varying trends, varying treatment-by-time interactions
Schutte_RML5 <- suppressWarnings(
  update(Schutte_RML4, random = ~ week + trt_week | case,
         control=lmeControl(msMaxIter = 50, apVar=FALSE, returnObject=TRUE))
)

# mod <- Laski_iid
# grps <- Laski$case
# invert <- TRUE
# sigma_scale <- TRUE
# R_list <- build_corr_mats(mod)

test_that("targetVariance() works with 2-level models.", {
  test_Sigma_mats(Laski_iid, Laski$case)
  test_Sigma_mats(Laski_het, Laski$case)
  test_Sigma_mats(Laski_AR1, Laski$case)
  test_Sigma_mats(Laski_hetAR1, Laski$case)
  test_Sigma_mats(Laski_CAR1, Laski$case)
  test_Sigma_mats(Laski_MA1, Laski$case)
  test_Sigma_mats(Laski_MA2, Laski$case)
  test_Sigma_mats(Laski_AR1MA1, Laski$case)
  test_Sigma_mats(Bodyweight_Gaus)
  test_Sigma_mats(Lambert_RML, subset(Lambert, measure=="academic response")$case)
  test_Sigma_mats(Anglesea_RML, Anglesea$case)
  test_Sigma_mats(Saddler_quality_RML, subset(Saddler, measure=="writing quality")$case )
  test_Sigma_mats(Schutte_RML3, Schutte$case)
  test_Sigma_mats(Schutte_RML4, Schutte$case)
  test_Sigma_mats(Schutte_RML5, Schutte$case)
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
  expect_error(test_deriv_dims(Bodyweight_Gaus))
  test_deriv_dims(Lambert_RML)
  test_deriv_dims(Anglesea_RML)
  test_deriv_dims(Saddler_quality_RML)
  test_deriv_dims(Schutte_RML3)
  test_deriv_dims(Schutte_RML4)
  test_deriv_dims(Schutte_RML5)
})

test_that("Information matrices work with FIML too.", {
  test_with_FIML(Laski_iid)
  test_with_FIML(Laski_het)
  test_with_FIML(Laski_AR1)
  test_with_FIML(Laski_hetAR1)
  test_with_FIML(Laski_CAR1)
  test_with_FIML(Laski_MA1)
  test_with_FIML(Lambert_RML)
  test_with_FIML(Anglesea_RML)
  test_with_FIML(Saddler_quality_RML)
  test_with_FIML(Schutte_RML3)
  test_with_FIML(Schutte_RML4)
  test_with_FIML(Schutte_RML5)
})

test_that("dR_dcorStruct.corCAR1 returns the same result as dR_dcorStruct.corAR1.", {
  expect_equal(dR_dcorStruct.corCAR1(Laski_CAR1$modelStruct$corStruct),
               dR_dcorStruct.corAR1(Laski_AR1$modelStruct$corStruct),
               tol = 10^-5)
})

test_that("lmeinfo::g_mlm returns the same result as scdhlm::g_REML.", {

  skip_if_not_installed("scdhlm", minimum_version = "0.4.2")

  check_against_scdhlm(Laski_AR1,
                       p_lmeInfo = c(0,1), r_lmeInfo = c(1,0,0,0,1),
                       p_scdhlm = c(0,1), r_scdhlm = c(1,0,1,0,0))
})


test_that("Results do not depend on order of data.", {

  skip_on_cran()

  test_after_shuffling(Laski_iid, seed = 20)
  test_after_shuffling(Laski_het, seed = 17) # 20
  test_after_shuffling(Laski_AR1, seed = 20)
  test_after_shuffling(Laski_hetAR1, test = "full", seed = 17) # 20
  test_after_shuffling(Laski_CAR1, seed = 20)
  test_after_shuffling(Laski_MA1, seed = 20)
  test_after_shuffling(Lambert_RML, seed = 20)
  test_after_shuffling(Anglesea_RML, seed = 20)
  test_after_shuffling(Saddler_quality_RML, seed = 20)
  test_after_shuffling(Schutte_RML3, seed = 20)
  test_after_shuffling(Schutte_RML4, seed = 20)
})


test_that("Info matrices work with dropped observations.", {

  skip_on_cran()

  test_after_deleting(Laski_iid)
  test_after_deleting(Laski_het)
  test_after_deleting(Laski_AR1)
  test_after_deleting(Laski_hetAR1)
  test_after_deleting(Laski_CAR1)
  test_after_deleting(Laski_MA1)
  test_after_deleting(Lambert_RML)
  test_after_deleting(Anglesea_RML)
  test_after_deleting(Saddler_quality_RML)
  test_after_deleting(Schutte_RML3)
  test_after_deleting(Schutte_RML4)
})

test_that("New REML calculations work.", {

  check_REML2(Laski_iid)
  check_REML2(Laski_het)
  check_REML2(Laski_AR1)
  check_REML2(Laski_hetAR1)
  check_REML2(Laski_CAR1)
  check_REML2(Laski_MA1)
  check_REML2(Lambert_RML)
  check_REML2(Anglesea_RML)
  check_REML2(Saddler_quality_RML)
  check_REML2(Schutte_RML3)
  check_REML2(Schutte_RML4)
  check_REML2(Schutte_RML5)

})
