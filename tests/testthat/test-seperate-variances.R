library(nlme)
data(Laski, package = "scdhlm")


Laski$trt_rev <- factor(ifelse(Laski$treatment == "treatment", "a_treatment","b_baseline"), levels = c("a_treatment","b_baseline"))

test_that("The separate_variances option works for gls() models.", {

  # gls
  Laski_AR1_gls <- gls(outcome ~ 0 + case + case:treatment,
                       correlation = corAR1(0.2, ~ time | case),
                       data = Laski)

  Laski_het_gls <- gls(outcome ~ 0 + case + case:treatment,
                       weights = varIdent(form = ~ 1 | treatment),
                       data = Laski)

  AR1_gls_no_sep <- extract_varcomp(Laski_AR1_gls, separate_variances = FALSE)
  expect_warning(AR1_gls_sep <- extract_varcomp(Laski_AR1_gls, separate_variances = TRUE))
  expect_equal(AR1_gls_no_sep, AR1_gls_sep)

  het_gls_no_sep <- extract_varcomp(Laski_het_gls, separate_variances = FALSE)
  het_gls_sep <- extract_varcomp(Laski_het_gls, separate_variances = TRUE)
  expect_equal(names(het_gls_no_sep), c("cor_params", "var_params", "sigma_sq"))
  expect_equal(names(het_gls_sep), c("cor_params", "sigma_sq"))
  expect_equal(het_gls_no_sep$cor_params, het_gls_sep$cor_params)

  vcov_gls_sep <- varcomp_vcov(Laski_het_gls, separate_variances = TRUE)
  vcov_gls_nosep <- varcomp_vcov(Laski_het_gls, separate_variances = FALSE)
  expect_equal(vcov_gls_sep[1,1], vcov_gls_nosep[2,2])

  # Laski_rev <- Laski[order(Laski$case, Laski$trt_rev),]
  #
  # Laski_rev_gls <- gls(outcome ~ 0 + case + case:trt_rev,
  #                      weights = varIdent(form = ~ 1 | trt_rev),
  #                      data = Laski_rev)
  #
  # hte_gls_no_sep_rev <- extract_varcomp(Laski_rev_gls, separate_variances = FALSE)
  # hte_gls_sep_rev <- extract_varcomp(Laski_rev_gls, separate_variances = TRUE)
  # expect_equivalent(rev(het_gls_sep$sigma_sq), hte_gls_sep_rev$sigma_sq)
  # expect_equal(hte_gls_no_sep_rev$sigma_sq, het_gls_sep$sigma_sq[["treatment"]])

  data(Orthodont)
  Ortho_power <- gls(distance ~ Subject:age + Sex,
                     weights = varPower(),
                     data = Orthodont)
  expect_warning(extract_varcomp(Ortho_power, separate_variances = TRUE))
  expect_warning(varcomp_vcov(Ortho_power, separate_variances = TRUE))

})


test_that("The separate_variances option works for two-level lme() models.", {

  Laski_rev <- Laski[order(Laski$case, Laski$trt_rev),]

  Laski_AR1_lme <- lme(fixed = outcome ~ treatment,
                       random = ~ treatment | case,
                       correlation = corAR1(0.2, ~ time | case),
                       data = Laski)

  Laski_het_lme <- lme(fixed = outcome ~ treatment,
                       random = ~ treatment | case,
                       weights = varIdent(form = ~ 1 | treatment),
                       data = Laski)

  Laski_rev_lme <- lme(fixed = outcome ~ trt_rev,
                       random = ~ trt_rev | case,
                       weights = varIdent(form = ~ 1 | trt_rev),
                       data = Laski_rev)

  AR1_lme_no_sep <- extract_varcomp(Laski_AR1_lme, separate_variances = FALSE)
  expect_warning(AR1_lme_sep <- extract_varcomp(Laski_AR1_lme, separate_variances = TRUE))
  expect_equal(AR1_lme_no_sep, AR1_lme_sep)

  het_lme_no_sep <- extract_varcomp(Laski_het_lme, separate_variances = FALSE)
  het_lme_sep <- extract_varcomp(Laski_het_lme, separate_variances = TRUE)

  expect_equal(names(het_lme_no_sep), c("Tau", "cor_params", "var_params", "sigma_sq"))
  expect_equal(names(het_lme_sep), c("Tau", "cor_params", "sigma_sq"))
  expect_equal(het_lme_no_sep$Tau, het_lme_sep$Tau)
  expect_equal(het_lme_no_sep$cor_params, het_lme_sep$cor_params)
  expect_equal(
    het_lme_no_sep$sigma_sq * c(1, het_lme_no_sep$var_params^2),
    as.numeric(het_lme_sep$sigma_sq)
  )

  hte_lme_no_sep_rev <- extract_varcomp(Laski_rev_lme, separate_variances = FALSE)
  hte_lme_sep_rev <- extract_varcomp(Laski_rev_lme, separate_variances = TRUE)
  expect_equivalent(rev(het_lme_sep$sigma_sq), hte_lme_sep_rev$sigma_sq, tolerance = 1e-5)
  expect_equal(hte_lme_no_sep_rev$sigma_sq, het_lme_sep$sigma_sq[["treatment"]], tolerance = 1e-5)

  vcov_lme_sep <- varcomp_vcov(Laski_het_lme, separate_variances = TRUE)
  vcov_lme_nosep <- varcomp_vcov(Laski_het_lme, separate_variances = FALSE)
  expect_equal(vcov_lme_sep["sigma_sq.baseline","sigma_sq.baseline"], vcov_lme_nosep["sigma_sq","sigma_sq"])

  vcov_lme_sep_rev <- varcomp_vcov(Laski_rev_lme, separate_variances = TRUE)
  vcov_lme_nosep_rev <- varcomp_vcov(Laski_rev_lme, separate_variances = FALSE)
  expect_equal(vcov_lme_sep["sigma_sq.treatment","sigma_sq.treatment"], vcov_lme_nosep_rev["sigma_sq","sigma_sq"], tolerance = 1e-5)
  expect_equivalent(vcov_lme_sep[4:5,4:5], vcov_lme_sep_rev[5:4,5:4], tolerance = 1e-5)

})

test_that("The separate_variances option works for three-level models with multiple level-1 variance estimates.", {

  data(Thiemann2001, package = "scdhlm")

  Thiemann <- lme(fixed = outcome ~ time_c + treatment + trt_time,
                  random = ~ 1 | case/series,
                  correlation = corAR1(0, ~ time_c | case/series),
                  weights = varIdent(form = ~ 1 | treatment),
                  data = Thiemann2001)

  # extract_varcomp()
  Thiemann_no_sep <- extract_varcomp(Thiemann, separate_variances = FALSE)
  Thiemann_sep <- extract_varcomp(Thiemann, separate_variances = TRUE)
  expect_equal(names(Thiemann_no_sep), c("Tau", "cor_params", "var_params", "sigma_sq"))
  expect_equal(names(Thiemann_sep), c("Tau", "cor_params", "sigma_sq"))

  # g_mlm()
  g_Thiemann <- g_mlm(Thiemann, p_const = c(0,0,1,22), r_const = c(1,1,0,0,1))
  g_Thiemann_bs <- g_mlm(Thiemann, p_const = c(0,0,1,22), r_const = c(1,1,0,1,0), separate_variances = TRUE)
  g_Thiemann_trt <- g_mlm(Thiemann, p_const = c(0,0,1,22), r_const = c(1,1,0,0,1), separate_variances = TRUE)

  expect_equal(as.numeric(g_Thiemann$SE_theta[5]), as.numeric(g_Thiemann_bs$SE_theta[4]))
  expect_equal(g_Thiemann_bs$SE_theta, g_Thiemann_trt$SE_theta)

  expect_equal(g_Thiemann$delta_AB, g_Thiemann_bs$delta_AB)
  expect_equal(g_Thiemann$g_AB, g_Thiemann_bs$g_AB)
  expect_equal(g_Thiemann$SE_g_AB, g_Thiemann_bs$SE_g_AB)
  expect_equal(g_Thiemann$nu, g_Thiemann_bs$nu)

  expect_true(g_Thiemann_trt$delta_AB < g_Thiemann_bs$delta_AB)
  expect_true(g_Thiemann_trt$g_AB < g_Thiemann_bs$g_AB)


})

