skip_if_not_installed("mlmRev")

library(nlme)

data(egsingle, package = "mlmRev")

# Take a subsample of students
set.seed(20221019)
childIDs <- unique(egsingle$childid)
child_subset <- sample(childIDs, size = 500L)
egsingle <- subset(egsingle, childid %in% child_subset)

# clean up predictors
egsingle$fem <- 1L * (egsingle$female == "Female")
egsingle$retained <- as.integer(egsingle$retained)
egsingle$black <- as.integer(egsingle$black)
egsingle$hispanic <- as.integer(egsingle$hispanic)

tol <- 1e-3

test_that("Separate school random effects equivalent to pdDiag specification", {
  skip_on_cran()

  m1 <- lme(math ~ year + size + female + black + hispanic,
            random = list(schoolid = pdDiag(~ 1 + year + fem),
                          ~ year | childid),
            data = egsingle, method = "ML")


  varcomp1 <- extract_varcomp(m1, vector = TRUE)
  vcvcov1 <- varcomp_vcov(m1)

  m2 <- lme(math ~ year + size + female + black + hispanic,
            random = list(~ 1 | schoolid, ~ 0 + year | schoolid, ~ 0 + fem | schoolid,
                          ~ year | childid),
            data = egsingle, method = "ML")

  expect_equivalent(
    varcomp1,
    extract_varcomp(m2, vector = TRUE),
    tol = tol
  )

  expect_equivalent(
    sqrt(diag(vcvcov1)),
    sqrt(diag(varcomp_vcov(m2))),
    tol = tol
  )

  m3 <- lme(math ~ year + size + female + black + hispanic,
            random = list(schoolid = pdDiag(~ 1 + year),
                          ~ 0 + fem | schoolid,
                          ~ year | childid),
            data = egsingle, method = "ML")

  expect_equivalent(
    varcomp1,
    extract_varcomp(m3, vector = TRUE),
    tol = tol
  )

  expect_equivalent(
    sqrt(diag(vcvcov1)),
    sqrt(diag(varcomp_vcov(m3))),
    tol = tol
  )

  m4 <- lme(math ~ year + size + female + black + hispanic,
            random = list(schoolid = pdDiag(~ 1 + fem),
                          schoolid = pdSymm(~ 0 + year),
                          ~ year | childid),
            data = egsingle, method = "ML")

  ord <- c(1,3,2,4:7)

  expect_equivalent(
    varcomp1,
    extract_varcomp(m4, vector = TRUE)[ord],
    tol = tol
  )

  expect_equivalent(
    sqrt(diag(vcvcov1)),
    sqrt(diag(varcomp_vcov(m4)))[ord],
    tol = tol
  )

})



test_that("Separate child random effects equivalent to pdDiag specification", {
  skip_on_cran()

  m1 <- lme(math ~ year + size + female + black + hispanic + retained,
            random = list(schoolid = pdDiag(~ 1 + year),
                          childid = pdDiag( ~ 1 + year + retained)),
            data = egsingle, method = "ML",
            control = lmeControl(tolerance = 1e-8))

  varcomp1 <- extract_varcomp(m1, vector = TRUE)
  vcvcov1 <- varcomp_vcov(m1)

  m2 <- lme(math ~ year + size + female + black + hispanic + retained,
            random = list(schoolid = pdDiag(~ 1 + year),
                          ~ 1 | childid, ~ 0 + year | childid, ~ 0 + retained | childid),
            data = egsingle, method = "ML",
            control = lmeControl(tolerance = 1e-8))

  expect_equivalent(
    varcomp1,
    extract_varcomp(m2, vector = TRUE),
    tol = tol
  )

  expect_equivalent(
    sqrt(diag(vcvcov1)),
    sqrt(diag(varcomp_vcov(m2))),
    tol = tol
  )

  m3 <- lme(math ~ year + size + female + black + hispanic + retained,
            random = list(schoolid = pdDiag(~ 1 + year),
                          childid = pdDiag(~ 0 + year + retained),
                          ~ 1 | childid),
            data = egsingle, method = "ML",
            control = lmeControl(tolerance = 1e-8))

  ord3 <- c(1,2,5,3,4,6)
  varcomp3 <- extract_varcomp(m3, vector = TRUE)

  expect_equivalent(
    varcomp1,
    varcomp3[ord3],
    tol = tol
  )

  expect_equivalent(
    sqrt(diag(vcvcov1)),
    sqrt(diag(varcomp_vcov(m3)))[ord3],
    tol = tol
  )


  m4 <- lme(math ~ year + size + female + black + hispanic + retained,
            random = list(schoolid = pdDiag(~ 1 + year),
                          childid = pdDiag(~ 1 + retained), ~ 0 + year | childid),
            data = egsingle, method = "ML",
            control = lmeControl(tolerance = 1e-8))

  ord4 <- c(1,2,3,5,4,6)

  expect_equivalent(
    varcomp3[ord3],
    extract_varcomp(m4, vector = TRUE)[ord4],
    tol = tol
  )

  expect_equivalent(
    sqrt(diag(vcvcov1)),
    sqrt(diag(varcomp_vcov(m4)))[ord4],
    tol = tol
  )

})


test_that("Separate order of same-level random effects does not matter.", {


  m1 <- lme(math ~ year + size + female + black + hispanic,
            random = list(schoolid = pdDiag(~ 1 + year),
                          schoolid = pdSymm(~ 0 + black + hispanic),
                          ~ year | childid),
            data = egsingle, method = "REML")

  varcomp1 <- extract_varcomp(m1, vector = TRUE)
  vcvcov1 <- varcomp_vcov(m1)

  m2 <- lme(math ~ year + size + female + black + hispanic,
            random = list(schoolid = pdSymm(~ 0 + black + hispanic),
                          schoolid = pdDiag(~ 1 + year),
                          ~ year | childid),
            data = egsingle, method = "REML")

  varcomp2 <- extract_varcomp(m2, vector = TRUE)
  vcvcov2 <- varcomp_vcov(m2)

  ord <- c(4,5,1:3,6:9)

  expect_equivalent(
    varcomp1,
    varcomp2[ord],
    tol = tol
  )

  expect_equivalent(
    sqrt(diag(vcvcov1)),
    sqrt(diag(vcvcov2))[ord],
    tol = tol
  )

})

