library(nlme)

# Thiemann 2001

data(Thiemann2001)

Thiemann2001_RML1 <- lme(fixed = outcome ~ treatment,
                     random = ~ 1 | case/series,
                     correlation = corAR1(0, ~ time | case/series),
                     data = Thiemann2001)

Thiemann2001_RML2 <- lme(fixed = outcome ~ treatment,
                      random = ~ treatment | case/series,
                      correlation = corAR1(0, ~ time | case/series),
                      data = Thiemann2001,
                      control=lmeControl(msMaxIter = 200, apVar=FALSE, returnObject=TRUE))

Thiemann2001_RML3 <- lme(fixed = outcome ~ time_c + treatment + trt_time,
                      random = ~ 1 | case/series,
                      correlation = corAR1(0, ~ time_c | case/series),
                      data = Thiemann2001)

Thiemann2001_RML4 <- lme(fixed = outcome ~ time_c + treatment + trt_time,
                      random = list(~ 1 | case, ~ treatment | series),
                      correlation = corAR1(0, ~ time_c | case/series),
                      data = Thiemann2001,
                      control=lmeControl(msMaxIter = 50, apVar=FALSE, returnObject=TRUE))


# Thiemann 2004

data(Thiemann2004)

Thiemann2004_RML1 <- lme(fixed = outcome ~ treatment,
                     random = ~ 1 | case/series,
                     correlation = corAR1(0, ~ time | case/series),
                     data = Thiemann2004)

Thiemann2004_RML2 <- lme(fixed = outcome ~ treatment,
                         random =list(~ 1 | case, ~ treatment | series),
                         correlation = corAR1(0, ~ time | case/series),
                         data = Thiemann2004,
                         control=lmeControl(msMaxIter = 50, apVar=FALSE, returnObject=TRUE))

Thiemann2004_RML3 <- lme(fixed = outcome ~ time_c + treatment + trt_time,
                      random = ~ 1 | case/series,
                      correlation = corAR1(0, ~ time_c | case/series),
                      data = Thiemann2004)

Thiemann2004_RML4 <- lme(fixed = outcome ~ time_c + treatment + trt_time,
                      random = list(~ treatment | case, ~ time_c | series),
                      correlation = corAR1(0, ~ time_c | case/series),
                      data = Thiemann2004,
                      control=lmeControl(msMaxIter = 50, apVar=FALSE, returnObject=TRUE))


# Bryant 2016

data(Bryant2016)

Bryant2016_RML1 <- lme(fixed = outcome ~ treatment,
                    random = ~ 1 | school/case,
                    correlation = corAR1(0, ~ session | school/case),
                    data = Bryant2016)

Bryant2016_RML2 <- lme(fixed = outcome ~ treatment,
                     random = ~ treatment | school/case,
                     correlation = corAR1(0, ~ session | school/case),
                     data = Bryant2016)

Bryant2016_RML3 <- lme(fixed = outcome ~ session_c + treatment + trt_time,
                     random = ~ 1 | school/case,
                     correlation = corAR1(0, ~ session_c | school/case),
                     data = Bryant2016)

Bryant2016_RML4 <- lme(fixed = outcome ~ session_c + treatment + trt_time,
                       random = list(~ 1 | school, ~ session_c | case),
                       correlation = corAR1(0, ~ session_c | school/case),
                       data = Bryant2016,
                       control=lmeControl(msMaxIter = 50, apVar=FALSE, returnObject=TRUE))

# Bryant 2018
data(Bryant2018)

Bryant2018_RML1 <- lme(fixed = outcome ~ treatment,
                     random = ~ 1 | school/case,
                     correlation = corAR1(0, ~ session | school/case),
                     data = Bryant2018)

Bryant2018_RML2 <- lme(fixed = outcome ~ treatment,
                     random = ~ treatment | school/case,
                     correlation = corAR1(0, ~ session | school/case),
                     data = Bryant2018)

Bryant2018_RML3 <- lme(fixed = outcome ~ session_c + treatment + session_trt,
                     random = ~ 1 | school/case,
                     correlation = corAR1(0, ~ session_c | school/case),
                     data = Bryant2018)

Bryant2018_RML4 <- lme(fixed = outcome ~ session_c + treatment + session_trt,
                     random = ~ session_c | school/case,
                     correlation = corAR1(0, ~ session_c | school/case),
                     data = Bryant2018,
                     control=lmeControl(msMaxIter = 50, apVar=FALSE, returnObject=TRUE))


test_that("targetVariance() works with 3-level models.", {

  test_Sigma_mats(Thiemann2001_RML1, Thiemann2001$case)
  test_Sigma_mats(Thiemann2001_RML2, Thiemann2001$case)
  test_Sigma_mats(Thiemann2001_RML3, Thiemann2001$case)
  test_Sigma_mats(Thiemann2001_RML4, Thiemann2001$case)
  test_Sigma_mats(Thiemann2004_RML1, Thiemann2004$case)
  test_Sigma_mats(Thiemann2004_RML2, Thiemann2004$case)
  test_Sigma_mats(Thiemann2004_RML3, Thiemann2004$case)
  test_Sigma_mats(Thiemann2004_RML4, Thiemann2004$case)
  test_Sigma_mats(Bryant2016_RML1, Bryant2016$school)
  test_Sigma_mats(Bryant2016_RML2, Bryant2016$school)
  test_Sigma_mats(Bryant2016_RML3, Bryant2016$school)
  test_Sigma_mats(Bryant2016_RML4, Bryant2016$school)
  test_Sigma_mats(Bryant2018_RML1, Bryant2018$school)
  test_Sigma_mats(Bryant2018_RML2, Bryant2018$school)
  test_Sigma_mats(Bryant2018_RML3, Bryant2018$school)
  test_Sigma_mats(Bryant2018_RML4, Bryant2018$school)

  expect_error(test_Sigma_mats(Bryant2016_RML1, Thiemann2001$case))
  expect_error(test_Sigma_mats(Bryant2016_RML1, Bryant2016$case))
  expect_error(test_Sigma_mats(Bryant2016_RML1, Bryant2018$school))

  expect_error(test_Sigma_mats(Bryant2018_RML1, Thiemann2001$case))
  expect_error(test_Sigma_mats(Bryant2018_RML1, Bryant2018$case))
  expect_error(test_Sigma_mats(Bryant2018_RML1, Bryant2016$school))

})
