context("Examples from scdhlm")

skip_if_not_installed("scdhlm")

library(scdhlm)

data(Lambert)
Lambert_RML <- lme(fixed = outcome ~ treatment,
                    random = ~ 1 | case,
                    correlation = corAR1(0.1, ~ time | case),
                    data = subset(Lambert, measure=="academic response"))


data(Anglesea)
Anglesea_RML <- lme(fixed = outcome ~ condition,
                    random = ~ condition | case,
                    correlation = corAR1(0.2, ~ session | case),
                    data = Anglesea)


data(Saddler)
Saddler_quality_RML <- lme(fixed = outcome ~ treatment,
                           random = ~ 1 | case,
                           correlation = corAR1(0, ~ time | case),
                           data = subset(Saddler, measure=="writing quality"))


data(Laski)
# create trt-by-time interaction
Laski$trt_time <- with(Laski,
                       unlist(tapply((treatment=="treatment") * time,
                       list(treatment, case),
                       function(x) x - min(x))) + (treatment=="treatment"))
# time-point constants
A_Laski <- 4
B_Laski <- 12
# center at follow-up time (time 12)
Center_Laski <- B_Laski
Laski$time <- Laski$time - Center_Laski

# Varying intercepts, fixed treatment effect, fixed trends
Laski_RML3 <- lme(fixed = outcome ~ time + treatment + trt_time,
                  random = ~ 1 | case,
                  correlation = corAR1(0.1, ~ time | case),
                  data = Laski)

# Varying intercepts, fixed treatment effect, varying trends

Laski_RML4 <- suppressWarnings(
  update(Laski_RML3, random = ~ time | case,
         control=lmeControl(msMaxIter = 50, apVar=FALSE, returnObject=TRUE))
)


data(Schutte)
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


test_that("lmeInfo::g_mlm returns the same result as scdhlm::g_REML.", {

  skip_if_not_installed("scdhlm", minimum_version = "0.4.2")

  check_against_scdhlm(Lambert_RML,
                       p_lmeInfo = c(0,1), r_lmeInfo = c(1,0,1), r_scdhlm = c(1,0,1))

  check_against_scdhlm(Anglesea_RML,
                       p_lmeInfo = c(0,1), r_lmeInfo = c(1,0,0,0,1), r_scdhlm = c(1,0,1,0,0))

  check_against_scdhlm(Saddler_quality_RML,
                       p_lmeInfo = c(0,1), r_lmeInfo = c(1,0,1), r_scdhlm = c(1,0,1))

  # Laski
  check_against_scdhlm(Laski_RML3,
                       p_lmeInfo = c(0,0,1,8), r_lmeInfo = c(1,0,1), r_scdhlm = c(1,0,1))

  check_against_scdhlm(Laski_RML4,
                       p_lmeInfo = c(0,0,1,8), r_lmeInfo = c(1,0,0,0,1), r_scdhlm = c(1,0,1,0,0))

  # Schutte
  check_against_scdhlm(Schutte_RML3,
                       p_lmeInfo = c(0,0,1,7), r_lmeInfo = c(1,0,1), r_scdhlm = c(1,0,1))

  check_against_scdhlm(Schutte_RML4,
                       p_lmeInfo = c(0,0,1,7), r_lmeInfo = c(1,0,0,0,1), r_scdhlm = c(1,0,1,0,0))

  check_against_scdhlm(Schutte_RML5,
                       p_lmeInfo = c(0,0,1,7),
                       r_lmeInfo = c(1,0,0,0,0,0,0,1), r_scdhlm = c(1,0,1,0,0,0,0,0))

})

test_that("lmeInfo::CI_g returns the correct CIs for Schutte examples.", {

  skip_if_not_installed("scdhlm", minimum_version = "0.4.2")

  g_RML4 <- g_mlm(Schutte_RML4, p_const = c(0,0,1,7), r_const = c(1,0,0,0,1))
  g_RML4_scdhlm <- suppressWarnings(scdhlm::g_REML(Schutte_RML4, p_const = c(0,0,1,7), r_const = c(1,0,1,0,0)))
  g_RML5 <- g_mlm(Schutte_RML5, p_const = c(0,0,1,7), r_const = c(1,0,0,0,0,0,0,1))
  g_RML5_scdhlm <- suppressWarnings(scdhlm::g_REML(Schutte_RML5, p_const = c(0,0,1,7), r_const = c(1,0,1,0,0,0,0,0)))

  # symmetric CIs

  CI_g_RML4 <- round(CI_g(g_RML4, symmetric = TRUE), 2)
  expect_equal(-2.02, CI_g_RML4[1])
  expect_equal(-.01, CI_g_RML4[2])

  CI_g_RML5 <- round(CI_g(g_RML5, symmetric = TRUE), 2)
  expect_equal(-3.29, CI_g_RML5[1])
  expect_equal(.64, CI_g_RML5[2])

  # approximate non-central t confidence interval
  appro_CI4_lmeInfo <- CI_g(g_RML4, symmetric = FALSE)
  appro_CI5_lmeInfo <- CI_g(g_RML5, symmetric = FALSE)

  appro_CI4_scdhlm <- suppressWarnings(scdhlm::CI_g(g_RML4_scdhlm, symmetric = FALSE))
  appro_CI5_scdhlm <- suppressWarnings(scdhlm::CI_g(g_RML5_scdhlm, symmetric = FALSE))

  expect_equal(appro_CI4_lmeInfo, appro_CI4_scdhlm)
  expect_equal(appro_CI5_lmeInfo, appro_CI5_scdhlm)

})

test_that("lmeInfo::summary() and lmeInfo::print() return output.", {

  skip_if_not_installed("scdhlm", minimum_version = "0.4.2")

  g_RML4 <- g_mlm(Schutte_RML4, p_const = c(0,0,1,7), r_const = c(1,0,0,0,1))
  expect_output(summary(g_RML4))
  expect_output(print(g_RML4))

  g_RML5 <- g_mlm(Schutte_RML5, p_const = c(0,0,1,7), r_const = c(1,0,0,0,0,0,0,1))
  expect_output(summary(g_RML5))
  expect_output(print(g_RML5))

})
