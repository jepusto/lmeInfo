skip_if_not_installed("lme4")
skip_if_not_installed("Matrix")
skip_if_not_installed("merDeriv")

library(lme4)
library(nlme)
suppressWarnings(library(merDeriv))

vcov_params <- function(x, sigma_sq) {
  d <- sqrt(2 * length(x) + 1/4) - 1/2
  res <- matrix(0, nrow = d, ncol = d)
  lower_tri <- lower.tri(res, diag = TRUE)
  res[lower_tri] <- x
  sigma_sq * (res %*% t(res))[lower_tri]
}

data(sleepstudy)

test_that("lmeInfo results comparable to merDeriv for sleepstudy data, estimated by FML.", {

  # sleepstudy data

  sleep_lme4_FML <- lmer(Reaction ~ Days + (Days|Subject),
                         data = sleepstudy,
                         REML = FALSE)
  sleep_nlme_FML <- lme(Reaction ~ Days,
                        random = ~ Days | Subject,
                        data = sleepstudy,
                        method = "ML")

  vcov_merDeriv_FML <- vcov(sleep_lme4_FML,
                            information = "expected",
                            full = TRUE, ranpar = "var")

  vcov_lmeInfo_FML <- varcomp_vcov(sleep_nlme_FML)

  # Full information maximum likelihood
  sigma_sq_lme4_FML <- getME(sleep_lme4_FML, "sigma")^2
  Tau_lme4_FML <- vcov_params(getME(sleep_lme4_FML, "theta"), sigma_sq = sigma_sq_lme4_FML)
  vc_FML <- extract_varcomp(sleep_nlme_FML)
  expect_equal(sigma_sq_lme4_FML, vc_FML$sigma_sq, tol = 10^-4)
  expect_equal(Tau_lme4_FML, vc_FML$Tau$Subject, check.attributes = FALSE, tol = 10^-4)

  expect_equal(as.matrix(vcov_merDeriv_FML[1:2,1:2]),
               vcov(sleep_nlme_FML),
               check.attributes = FALSE, tol = 10^-4)
  expect_equal(as.matrix(vcov_merDeriv_FML[3:6,3:6]),
               vcov_lmeInfo_FML,
               check.attributes = FALSE, tol = 10^-4)

})



test_that("lmeInfo results comparable to merDeriv for sleepstudy data, estimated by REML.", {

  # sleepstudy data

  sleep_lme4_REML <- lmer(Reaction ~ Days + (Days|Subject),
                          data = sleepstudy,
                          REML = TRUE)
  sleep_nlme_REML <- lme(Reaction ~ Days,
                         random = ~ Days | Subject,
                         data = sleepstudy,
                         method = "REML")

  vcov_merDeriv_REML <- vcov(sleep_lme4_REML,
                             information = "expected",
                             full = TRUE, ranpar = "var")

  vcov_lmeInfo_REML <- varcomp_vcov(sleep_nlme_REML)

  sigma_sq_lme4_REML <- getME(sleep_lme4_REML, "sigma")^2
  Tau_lme4_REML <- vcov_params(getME(sleep_lme4_REML, "theta"), sigma_sq = sigma_sq_lme4_REML)
  vc_REML <- extract_varcomp(sleep_nlme_REML)
  expect_equal(sigma_sq_lme4_REML, vc_REML$sigma_sq, tol = 10^-4)
  expect_equal(Tau_lme4_REML, vc_REML$Tau$Subject, check.attributes = FALSE, tol = 10^-3)

  expect_equal(as.matrix(vcov_merDeriv_REML[1:2,1:2]),
               vcov(sleep_nlme_REML),
               check.attributes = FALSE, tol = 10^-3)
  expect_equal(as.matrix(vcov_merDeriv_REML[3:6,3:6]),
               vcov_lmeInfo_REML,
               check.attributes = FALSE, tol = 10^-3)

})

data(bdf, package = "mlmRev")

test_that("lmeInfo results comparable to merDeriv for bdf data, estimated by FML.", {
  skip("Can't get the model parameter estimates to agree.")

  # bdf data

  bdf_lme4_FML <- suppressWarnings(
    lmer(langPOST ~ sex + Minority + aritPRET + (sex + aritPRET | schoolNR),
         data = bdf, REML = FALSE)
  )

  bdf_nlme_FML <- lme(langPOST ~ sex + Minority + aritPRET,
                      random = ~ sex + aritPRET | schoolNR,
                      data = bdf,
                      method = "ML")

  vcov_merDeriv_FML <- vcov(bdf_lme4_FML,
                            information = "expected",
                            full = TRUE, ranpar = "var")

  vcov_lmeInfo_FML <- varcomp_vcov(bdf_nlme_FML)

  sigma_sq_lme4_FML <- getME(bdf_lme4_FML, "sigma")^2
  Tau_lme4_FML <- vcov_params(getME(bdf_lme4_FML, "theta"), sigma_sq = sigma_sq_lme4_FML)
  vc_FML <- extract_varcomp(bdf_nlme_FML)
  expect_equal(sigma_sq_lme4_FML, vc_FML$sigma_sq, tol = 10^-3)
  expect_equal(Tau_lme4_FML[c(1,2,4,3,5,6)], vc_FML$Tau$schoolNR,
               check.attributes = FALSE, tol = 10^-2)

  expect_equal(as.matrix(vcov_merDeriv_FML[1:4,1:4]),
               vcov(bdf_nlme_FML),
               check.attributes = FALSE, tol = 5 * 10^-3)
  # This did not agree probably due to the order of parameters.
  expect_equal(as.matrix(vcov_merDeriv_FML[c(5,6,8,7,9:11),c(5,6,8,7,9:11)]),
               vcov_lmeInfo_FML,
               check.attributes = FALSE, tol = 10^-2)


})

test_that("lmeInfo results comparable to merDeriv for bdf data, estimated by REML.", {
  skip("Can't get the model parameter estimates to agree")

  bdf_lme4_REML <- suppressWarnings(
    lmer(langPOST ~ sex + Minority + aritPRET + (sex + aritPRET | schoolNR),
                        data = bdf, REML = TRUE))

  bdf_nlme_REML <- lme(langPOST ~ sex + Minority + aritPRET,
                       random = ~ sex + aritPRET | schoolNR,
                       data = bdf,
                       method = "REML")

  vcov_merDeriv_REML <- vcov(bdf_lme4_REML,
                             information = "expected",
                             full = TRUE, ranpar = "var")

  vcov_lmeInfo_REML <- varcomp_vcov(bdf_nlme_REML)

  sigma_sq_lme4_REML <- getME(bdf_lme4_REML, "sigma")^2
  Tau_lme4_REML <- vcov_params(getME(bdf_lme4_REML, "theta"), sigma_sq = sigma_sq_lme4_REML)
  vc_REML <- extract_varcomp(bdf_nlme_REML)
  expect_equal(sigma_sq_lme4_REML, vc_REML$sigma_sq, tol = 10^-3)
  expect_equal(Tau_lme4_REML, vc_REML$Tau$schoolNR,
               check.attributes = FALSE, tol = 10^-4)

  expect_equal(as.matrix(vcov_merDeriv_REML[1:4,1:4]),
               vcov(bdf_nlme_REML),
               check.attributes = FALSE, tol = 10^-4)
  expect_equal(as.matrix(vcov_merDeriv_REML[5:8,5:8]),
               vcov_lmeInfo_REML,
               check.attributes = FALSE, tol = 10^-4)

})
