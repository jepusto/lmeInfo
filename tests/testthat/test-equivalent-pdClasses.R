skip_if_not_installed("mlmRev")
skip_on_cran()

library(nlme)

data(egsingle, package = "mlmRev")
egsingle$fem <- 1L * (egsingle$female == "Female")


test_that("Separate school random effects equivalent to pdDiag specification", {

  m1 <- lme(math ~ year * size + female + black + hispanic,
            random = list(schoolid = pdDiag(~ 1 + year + fem),
                          ~ year | childid),
            data = egsingle, method = "ML")


  varcomp1 <- extract_varcomp(m1, vector = TRUE)
  vcvcov1 <- varcomp_vcov(m1)

  m2 <- lme(math ~ year * size + female + black + hispanic,
            random = list(~ 1 | schoolid, ~ 0 + year | schoolid, ~ 0 + fem | schoolid,
                          ~ year | childid),
            data = egsingle, method = "ML")

  expect_equivalent(
    varcomp1,
    extract_varcomp(m2, vector = TRUE),
    tol = 1e-4
  )

  expect_equivalent(
    vcvcov1,
    varcomp_vcov(m2),
    tol = 1e-4
  )

  m3 <- lme(math ~ year * size + female + black + hispanic,
            random = list(schoolid = pdDiag(~ 1 + year),
                          ~ 0 + fem | schoolid,
                          ~ year | childid),
            data = egsingle, method = "ML")

  expect_equivalent(
    varcomp1,
    extract_varcomp(m3, vector = TRUE),
    tol = 1e-4
  )

  expect_equivalent(
    vcvcov1,
    varcomp_vcov(m3),
    tol = 1e-4
  )

  m4 <- lme(math ~ year * size + female + black + hispanic,
            random = list(schoolid = pdDiag(~ 1 + fem),
                          schoolid = pdSymm(~ 0 + year),
                          ~ year | childid),
            data = egsingle, method = "ML")

  ord <- c(1,3,2,4:7)

  expect_equivalent(
    varcomp1,
    extract_varcomp(m4, vector = TRUE)[ord],
    tol = 1e-4
  )

  expect_equivalent(
    vcvcov1,
    varcomp_vcov(m4)[ord,ord],
    tol = 1e-4
  )
})



test_that("Separate child random effects equivalent to pdDiag specification", {

  m1 <- lme(math ~ year * size + female + black + hispanic,
            random = list(schoolid = pdDiag(~ 1 + year),
                          childid = pdDiag( ~ 1 + year + retained)),
            data = egsingle, method = "ML")

  varcomp1 <- extract_varcomp(m1, vector = TRUE)
  vcvcov1 <- varcomp_vcov(m1)

  m2 <- lme(math ~ year * size + female + black + hispanic,
            random = list(schoolid = pdDiag(~ 1 + year),
                          ~ 1 | childid, ~ 0 + year | childid, ~ 0 + retained | childid),
            data = egsingle, method = "ML")

  expect_equivalent(
    varcomp1,
    extract_varcomp(m2, vector = TRUE),
    tol = 1e-4
  )

  expect_equivalent(
    vcvcov1,
    varcomp_vcov(m2),
    tol = 1e-4
  )

  m3 <- lme(math ~ year * size + female + black + hispanic,
            random = list(schoolid = pdDiag(~ 1 + year),
                          ~ 1 | childid, ~ 0 + year | childid, ~ 0 + retained | childid),
            data = egsingle, method = "ML")

  expect_equivalent(
    varcomp1,
    extract_varcomp(m3, vector = TRUE),
    tol = 1e-4
  )

  expect_equivalent(
    vcvcov1,
    varcomp_vcov(m3),
    tol = 1e-4
  )

  m4 <- lme(math ~ year * size + female + black + hispanic,
            random = list(schoolid = pdDiag(~ 1 + year),
                          childid = pdDiag(~ 1 + retained), ~ 0 + year | childid),
            data = egsingle, method = "ML")

  ord <- c(1,3,2,4:7)

  expect_equivalent(
    varcomp1,
    extract_varcomp(m4, vector = TRUE)[ord],
    tol = 1e-4
  )

  expect_equivalent(
    vcvcov1,
    varcomp_vcov(m4)[ord,ord],
    tol = 1e-4
  )

})


test_that("Separate school random effects equivalent to pdDiag specification", {



})

