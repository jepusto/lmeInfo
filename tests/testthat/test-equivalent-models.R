library(nlme)
data(Thiemann2001)

Thiemann2001 <- Thiemann2001[sample(nrow(Thiemann2001)),]

Thiemann2001_3level <- lme(fixed = outcome ~ treatment,
                           random = ~ 1 | case/series,
                           data = Thiemann2001,
                           control = lmeControl(tolerance = 10^-8))

Thiemann2001_CScorr <- lme(fixed = outcome ~ treatment,
                           random = ~ 1 | case,
                           correlation = corCompSymm(form = ~ 1 | case / series),
                           data = Thiemann2001,
                           control = lmeControl(tolerance = 10^-8))

test_that("vcov matrices are equivalent for equivalent models.", {

  expect_equal(fixef(Thiemann2001_3level), fixef(Thiemann2001_CScorr), tol = 10^-6)

  vc_3level <- extract_varcomp(Thiemann2001_3level)
  vc_CScorr <- extract_varcomp(Thiemann2001_CScorr)
  expect_equal(vc_3level$Tau$case, vc_CScorr$Tau$case, tol = 10^-3)
  expect_equal(as.numeric(vc_3level$Tau$series) + vc_3level$sigma_sq, vc_CScorr$sigma_sq, tol = 10^-4)

  ICC_3level <- as.numeric(vc_3level$Tau$series / (vc_3level$Tau$series + vc_3level$sigma_sq))
  expect_equal(ICC_3level, vc_CScorr$cor_params, tol = 10^-3)

  vcov_3level_exp <- varcomp_vcov(Thiemann2001_3level)
  vcov_CScorr_exp <- varcomp_vcov(Thiemann2001_CScorr)

  expect_equal(vcov_3level_exp[1,1] / vcov_CScorr_exp[1,1], 1, tol = 10^-3)
  expect_equal(sum(vcov_3level_exp[2:3,2:3]) / vcov_CScorr_exp[3,3], 1, tol = 10^-3)
  expect_equal(sum(vcov_3level_exp[2:3,1]) / vcov_CScorr_exp[1,3], 1, tol = 10^-3)

})
