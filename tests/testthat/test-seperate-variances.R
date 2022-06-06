library(nlme)
data(Laski, package = "scdhlm")

Laski$trt_rev <- ifelse(Laski$treatment == "treatment", "a_treatment","b_baseline")


test_that("The separate_variances option works for gls() models.", {

  # gls

  Laski_AR1_gls <- gls(outcome ~ 0 + case + case:treatment,
                       correlation = corAR1(0.2, ~ time | case),
                       data = Laski)

  Laski_het_gls <- gls(outcome ~ 0 + case + case:treatment,
                       weights = varIdent(form = ~ 1 | treatment),
                       data = Laski)
  Laski_rev_gls <- gls(outcome ~ 0 + case + case:trt_rev,
                       weights = varIdent(form = ~ 1 | trt_rev),
                       data = Laski[order(Laski$case, Laski$trt_rev),])

  AR1_gls_no_sep <- extract_varcomp(Laski_AR1_gls, separate_variances = FALSE)
  expect_warning(AR1_gls_sep <- extract_varcomp(Laski_AR1_gls, separate_variances = TRUE))
  expect_equal(AR1_gls_no_sep, AR1_gls_sep)

  het_gls_no_sep <- extract_varcomp(Laski_het_gls, separate_variances = FALSE)
  het_gls_sep <- extract_varcomp(Laski_het_gls, separate_variances = TRUE)
  expect_equal(names(het_gls_no_sep), c("cor_params", "var_params", "sigma_sq"))
  expect_equal(names(het_gls_sep), c("cor_params", "sigma_sq"))
  expect_equal(het_gls_no_sep$cor_params, het_gls_sep$cor_params)

  data(Orthodont)
  Ortho_power <- gls(distance ~ Subject:age + Sex,
                     weights = varPower(),
                     data = Orthodont)
  expect_warning(extract_varcomp(Ortho_power, separate_variances = TRUE))

})


test_that("The separate_variances option works for two-level lme() models.", {

  # lme

  Laski_AR1_lme <- lme(fixed = outcome ~ treatment,
                       random = ~ treatment | case,
                       correlation = corAR1(0.2, ~ time | case),
                       data = Laski)

  Laski_het_lme <- lme(fixed = outcome ~ treatment,
                       random = ~ treatment | case,
                       weights = varIdent(form = ~ 1 | treatment),
                       data = Laski)

  Laski_rev_lme <- lme(fixed = outcome ~ treatment,
                       random = ~ treatment | case,
                       weights = varIdent(form = ~ 1 | trt_rev),
                       data = Laski)

  AR1_lme_no_sep <- extract_varcomp(Laski_AR1_lme, separate_variances = FALSE)
  expect_warning(AR1_lme_sep <- extract_varcomp(Laski_AR1_lme, separate_variances = TRUE))
  expect_equal(AR1_lme_no_sep, AR1_lme_sep)

  het_lme_no_sep <- extract_varcomp(Laski_het_lme, separate_variances = FALSE)
  het_lme_sep <- extract_varcomp(Laski_het_lme, separate_variances = TRUE)
  het_lme_rev <- extract_varcomp(Laski_rev_lme, separate_variances = FALSE)

  expect_equal(names(het_lme_no_sep), c("Tau", "cor_params", "var_params", "sigma_sq"))
  expect_equal(names(het_lme_sep), c("Tau", "cor_params", "sigma_sq"))
  expect_equal(het_lme_no_sep$Tau, het_lme_sep$Tau)
  expect_equal(het_lme_no_sep$cor_params, het_lme_sep$cor_params)
  expect_equal(
    het_lme_no_sep$sigma_sq * c(1, het_lme_no_sep$var_params^2),
    as.numeric(het_lme_sep$sigma_sq)
  )


})

test_that("The separate_variances option works for three-level models with multiple level-1 variance estimates.", {

  data(Thiemann2001, package = "scdhlm")
  Thiemann <- lme(fixed = outcome ~ time_c + treatment + trt_time,
                  random = ~ 1 | case/series,
                  correlation = corAR1(0, ~ time_c | case/series),
                  weights = varIdent(form = ~ 1 | treatment),
                  data = Thiemann2001)

  Thiemann_no_sep <- extract_varcomp(Thiemann, separate_variances = FALSE)
  Thiemann_sep <- extract_varcomp(Thiemann, separate_variances = TRUE)
  expect_equal(names(Thiemann_no_sep), c("Tau", "cor_params", "var_params", "sigma_sq"))
  expect_equal(names(Thiemann_sep), c("Tau", "cor_params", "sigma_sq"))

})
