skip_if_not_installed("nlme")
skip_if_not_installed("scdhlm")

library(nlme)
library(scdhlm)

test_that("The degrees of freedom are not sensitive to the choice of centering value for two-level models", {

  data(Laski)

  A <- 4
  B <- 13

  # create trt-by-time interaction
  Laski$trt_time <- with(Laski,
                         unlist(tapply((treatment=="treatment") * time,
                                       list(treatment, case),
                                       function(x) x - min(x))) + (treatment=="treatment"))
  # center at follow-up time (B)
  C1 <- B
  Laski$time1 <- Laski$time - C1

  # Varying intercepts, fixed treatment effect, varying trends
  Laski_RML_C1 <- suppressWarnings(lme(fixed = outcome ~ time1 + treatment + trt_time,
                                       random = ~ time1 | case,
                                       correlation = corAR1(0.01, ~ time1 | case),
                                       data = Laski,
                                       control = lmeControl(msMaxIter = 50, apVar = FALSE, returnObject = TRUE)))

  Laski_g_C1 <- suppressWarnings(g_mlm(Laski_RML_C1,
                                       p_const = c(0, 0, 1, B-A),
                                       r_const = c(1, 2*(B-C1), (B-C1)^2, 0, 1),
                                       returnModel = TRUE))

  # center at start of series
  C2 <- 0
  Laski$time2 <- Laski$time - C2

  Laski_RML_C2 <- suppressWarnings(lme(fixed = outcome ~ time2 + treatment + trt_time,
                                       random = ~ time2 | case,
                                       correlation = corAR1(0.01, ~ time2 | case),
                                       data = Laski,
                                       control = lmeControl(msMaxIter = 50, apVar = FALSE, returnObject = TRUE)))
  Laski_g_C2 <- suppressWarnings(g_mlm(Laski_RML_C2,
                                       p_const = c(0, 0, 1, B-A),
                                       r_const = c(1, 2*(B-C2), (B-C2)^2, 0, 1),
                                       returnModel = TRUE))

  expect_equal(Laski_g_C1$delta_AB, Laski_g_C2$delta_AB, tol = .005)
  expect_equal(Laski_g_C1$g_AB, Laski_g_C2$g_AB, tol = .005)
  expect_equal(Laski_g_C1$SE_g_AB, Laski_g_C2$SE_g_AB, tol = .001)
  expect_equal(Laski_g_C1$nu, Laski_g_C2$nu, tol = .005)

})

test_that("The degrees of freedom are not sensitive to the choice of centering value for Thiemann 2001 data.", {

  data(Thiemann2001)

  A <- 9
  B <- 15

  # center at B
  C1 <- B
  Thiemann2001$time1 <- Thiemann2001$time - C1

  Thi_RML_C1 <- suppressWarnings(lme(fixed = outcome ~ time1 + treatment + trt_time,
                                     random = list(case = ~ 1, series = ~ time1 + trt_time),
                                     correlation = corAR1(0.01, ~ time1 | case/series),
                                     data = Thiemann2001,
                                     control = lmeControl(msMaxIter = 50, apVar = FALSE, returnObject = TRUE)))
  Thi_g_C1 <- g_mlm(Thi_RML_C1,
                    p_const = c(0,0,1,B-A),
                    r_const = c(1,0,0,0,0,0,1,0,1),
                    returnModel = TRUE)

  # center at start of series
  C2 <- 0
  Thiemann2001$time2 <- Thiemann2001$time - C2

  Thi_RML_C2 <- suppressWarnings(lme(fixed = outcome ~ time2 + treatment + trt_time,
                                     random = list(case = ~ 1, series = ~ time2 + trt_time),
                                     correlation = corAR1(0.01, ~ time2 | case/series),
                                     data = Thiemann2001,
                                     control = lmeControl(msMaxIter = 50, apVar = FALSE, returnObject = TRUE)))
  Thi_g_C2 <- g_mlm(Thi_RML_C2,
                    p_const = c(0,0,1,B-A),
                    r_const = c(1,2*(B-C2),(B-C2)^2,0,0,0,1,0,1),
                    returnModel = TRUE)

  # check whether BCSMD estimates, SE, and dfs match for different centering values
  expect_equal(Thi_g_C1$delta_AB, Thi_g_C2$delta_AB, tol = .001)
  expect_equal(Thi_g_C1$g_AB, Thi_g_C2$g_AB, tol = .001)
  expect_equal(Thi_g_C1$SE_g_AB, Thi_g_C2$SE_g_AB, tol = .001)
  expect_equal(Thi_g_C1$nu, Thi_g_C2$nu, tol = .1)


})

test_that("The degrees of freedom are not sensitive to the choice of centering value for Bryant 2018 data.", {

  skip_if(packageVersion('scdhlm') <= '0.6.0')

  data(Bryant2018)

  A <- 5
  B <- 49

  # center at follow-up time B
  C1 <- B
  Bryant2018$session1 <-  Bryant2018$session - C1

  Bry_RML_C1 <- suppressWarnings(lme(fixed = outcome ~ session1 + treatment + session_trt,
                                     random = list(group = ~ 1, case = ~ session1 + session_trt),
                                     correlation = corAR1(0.01, ~ session1 | group/case),
                                     data = Bryant2018,
                                     control = lmeControl(msMaxIter = 50, apVar = FALSE, returnObject = TRUE)))
  Bry_g_C1 <- g_mlm(Bry_RML_C1,
                    p_const = c(0,0,1,B-A),
                    r_const = c(1,0,0,0,0,0,1,0,1),
                    returnModel = TRUE)

  # center at start of series
  C2 <- 0
  Bryant2018$session2 <-  Bryant2018$session - C2

  Bry_RML_C2 <- suppressWarnings(lme(fixed = outcome ~ session2 + treatment + session_trt,
                                     random = list(group = ~ 1, case = ~ session2 + session_trt),
                                     correlation = corAR1(0.01, ~ session2 | group/case),
                                     data = Bryant2018,
                                     control = lmeControl(msMaxIter = 50, apVar = FALSE, returnObject = TRUE)))
  Bry_g_C2 <- g_mlm(Bry_RML_C2,
                    p_const = c(0,0,1,B-A),
                    r_const = c(1,2*(B-C2),(B-C2)^2,0,0,0,1,0,1),
                    returnModel = TRUE)

  # check whether BCSMD estimates, SE, and dfs match for different centering values
  expect_equal(Bry_g_C1$delta_AB, Bry_g_C2$delta_AB, tol = .005)
  expect_equal(Bry_g_C1$g_AB, Bry_g_C2$g_AB, tol = .005)
  expect_equal(Bry_g_C1$SE_g_AB, Bry_g_C2$SE_g_AB, tol = .001)
  expect_equal(Bry_g_C1$nu, Bry_g_C2$nu, tol = .01)

})
