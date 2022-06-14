skip_if_not_installed("mlmRev")

library(nlme)

data(bdf, package = "mlmRev")

bdf_long <-
  bdf |>
  subset(
    schoolNR %in% levels(schoolNR)[1:12],
    select = c(schoolNR, pupilNR, sex, Minority, IQ.verb, IQ.perf, aritPRET)
  ) |>
  droplevels() |>
  reshape(
    direction = "long",
    idvar = "pupilNR",
    varying = c("IQ.verb","IQ.perf","aritPRET"),
    v.names = "score",
    timevar = "measure",
    times = c("IQ.verb","IQ.perf","aritPRET")
  )

bdf_long_order <- with(bdf_long, order(schoolNR, pupilNR, measure))
bdf_long <- bdf_long[bdf_long_order,]

bdf_MV2L <- lme(score ~ 0 + measure,
                random = ~ 1| schoolNR,
                corr = corCompSymm(form = ~ 1 | schoolNR / pupilNR),
                weights = varIdent(form = ~ 1 | measure),
                data = bdf_long,
                control = lmeControl(maxIter = 200, msMaxIter = 200, tolerance = 1e-8, opt = "optim", niterEM = 50))

bdf_MV3L <- lme(score ~ 0 + measure,
                random = ~ 1| schoolNR / pupilNR,
                corr = corCompSymm(form = ~ 1 | schoolNR / pupilNR),
                weights = varIdent(form = ~ 1 | measure),
                data = bdf_long,
                control = lmeControl(maxIter = 200, msMaxIter = 200, tolerance = 1e-8, opt = "optim", niterEM = 50))

# introduce random missing to bdf_long
set.seed(20200311)
bdf_long$school_id <- bdf_long$schoolNR
levels(bdf_long$school_id) <- LETTERS[1:nlevels(bdf_long$school_id)]

bdf_long_wm <-
  bdf_long |>
  within({
    row_index = rbinom(n = nrow(bdf_long), size = 1, prob = 0.9)
    score = ifelse(row_index == 1, score, as.numeric(NA))
    measure_id = as.integer(factor(measure))
  }) |>
  subset(
    !is.na(score),
    select = -row_index
  )

row_numbers <- 1:nrow(bdf_long_wm)

bdf_long_shuff <-
  bdf_long_wm |>
  within({
    index = ifelse(schoolNR == 1, row_numbers, sample(row_numbers))
    id = ifelse(schoolNR == 1, "A","B")
  })

bdf_shuff_order <- with(bdf_long_shuff, order(id, index))
bdf_long_shuff <- bdf_long_shuff[bdf_shuff_order,]

bdf_wm <- lme(score ~ 0 + measure,
              random = ~ 1| schoolNR / pupilNR,
              corr = corSymm(form = ~ measure_id | schoolNR / pupilNR),
              weights = varIdent(form = ~ 1 | measure_id),
              data = bdf_long_wm,
              control=lmeControl(msMaxIter = 100, apVar = FALSE, returnObject = TRUE))

bdf_2L_wm_shuff <- lme(score ~ 0 + measure,
                       random = ~ 1| schoolNR,
                       corr = corSymm(form = ~ measure_id | schoolNR / pupilNR),
                       weights = varIdent(form = ~ 1 | measure_id),
                       data = bdf_long_shuff,
                       control=lmeControl(msMaxIter = 100, apVar = FALSE, returnObject = TRUE))

bdf_3L_wm_shuff <- lme(score ~ 0 + measure,
                    random = ~ 1| schoolNR / pupilNR,
                    corr = corSymm(form = ~ measure_id | schoolNR / pupilNR),
                    weights = varIdent(form = ~ 1 | measure_id),
                    data = bdf_long_shuff,
                    control=lmeControl(msMaxIter = 100, apVar = FALSE, returnObject = TRUE))

bdf_wm_id <- lme(score ~ 0 + measure,
                  random = ~ 1| school_id / pupilNR,
                  corr = corSymm(form = ~ measure_id | school_id / pupilNR),
                  weights = varIdent(form = ~ 1 | measure_id),
                  data = bdf_long_shuff,
                  control=lmeControl(msMaxIter = 100, apVar = FALSE, returnObject = TRUE))

mod <- bdf_MV2L


test_that("targetVariance() works with multivariate models.", {

  skip_on_cran()

  test_Sigma_mats(bdf_MV2L, bdf_long$schoolNR)
  test_Sigma_mats(bdf_MV3L, bdf_long$schoolNR)
  test_Sigma_mats(bdf_wm, bdf_long_wm$schoolNR)
  test_Sigma_mats(bdf_wm_id, bdf_long_shuff$school_id)
  test_Sigma_mats(bdf_2L_wm_shuff, bdf_long_shuff$school_id)
  test_Sigma_mats(bdf_3L_wm_shuff, bdf_long_shuff$school_id)
})

test_that("Derivative matrices are of correct dimension with multivariate models.", {

  skip_on_cran()

  test_deriv_dims(bdf_MV2L)
  test_deriv_dims(bdf_MV3L)
  test_deriv_dims(bdf_wm)
  test_deriv_dims(bdf_wm_id)
  test_deriv_dims(bdf_2L_wm_shuff)
  test_deriv_dims(bdf_3L_wm_shuff)
})

test_that("Information matrices work with FIML too.", {

  skip_on_cran()

  test_with_FIML(bdf_MV2L)
  test_with_FIML(bdf_MV3L)
  test_with_FIML(bdf_wm)
  test_with_FIML(bdf_wm_id)
  test_with_FIML(bdf_2L_wm_shuff)
  test_with_FIML(bdf_3L_wm_shuff)
})

test_that("New REML calculations work.", {

  skip_on_cran()

  check_REML2(bdf_MV2L)
  check_REML2(bdf_MV3L)
  check_REML2(bdf_wm)
  check_REML2(bdf_wm_id)
  check_REML2(bdf_2L_wm_shuff)
  check_REML2(bdf_3L_wm_shuff)

})

test_that("Info matrices work with dropped observations.", {

  skip_on_cran()

  test_after_deleting(bdf_MV2L, seed = 10)
  test_after_deleting(bdf_MV3L, seed = 10)
  test_after_deleting(bdf_wm, seed = 20)
  test_after_deleting(bdf_wm_id, seed = 30)
  test_after_deleting(bdf_2L_wm_shuff, seed = 40)
  test_after_deleting(bdf_3L_wm_shuff, seed = 40)
})

test_that("Results do not depend on order of data.", {

  skip_on_cran()

  test_after_shuffling(bdf_MV2L, keep_rows = 72, seed = 22)
  test_after_shuffling(bdf_MV3L, keep_rows = 72, seed = 20)

  test_after_shuffling(bdf_wm, keep_rows = 72, seed = 17)
  test_after_shuffling(bdf_wm_id, keep_rows = 72, seed = 20)
  test_after_shuffling(bdf_2L_wm_shuff, keep_rows = 72, seed = 24)
  test_after_shuffling(bdf_3L_wm_shuff, keep_rows = 72, seed = 29)
})

test_that("The separate_variances option works with multivariate models.", {

  bdf_MVML_no_sep <- extract_varcomp(bdf_MV3L, separate_variances = FALSE)
  bdf_MVML_sep <- extract_varcomp(bdf_MV3L, separate_variances = TRUE)
  expect_equal(names(bdf_MVML_no_sep), c("Tau", "cor_params", "var_params", "sigma_sq"))
  expect_equal(names(bdf_MVML_sep), c("Tau", "cor_params", "sigma_sq"))

  expect_equal(names(bdf_MVML_sep$sigma_sq), unique(bdf_long$measure))

})

test_that("The Fisher_info() works correctly for multivariate models with a `varIdent()` variance structure.", {

  info_sep <- Fisher_info(bdf_MV3L, separate_variances = TRUE)

  # hand calculation by index
  info <- Fisher_info(bdf_MV3L, separate_variances = FALSE)
  theta <- extract_varcomp(bdf_MV3L, separate_variances = FALSE)
  theta_names <- vapply(strsplit(names(unlist(theta)), split = "[.]"),
                        function(x) paste(unique(x), collapse = "."), character(1L))
  theta_reparam <- extract_varcomp(bdf_MV3L, separate_variances = TRUE)
  theta_reparam_names <- vapply(strsplit(names(unlist(theta_reparam)), split = "[.]"),
                                function(x) paste(unique(x), collapse = "."), character(1L))
  r12 <- length(unlist(theta[c(1,2)]))
  r34 <- length(unlist(theta[c(3,4)]))
  r <- length(unlist(theta))
  Jac_1 <- Jac_inv_1 <- diag(1, r12)
  Jac_2 <- Jac_inv_2 <- matrix(0, nrow = r12, ncol = r34)
  Jac_3 <- Jac_inv_3 <- t(Jac_2)
  Jac_41 <- rep(0, length(unlist(theta[3])))
  Jac_42 <- 1
  Jac_43 <- diag(as.numeric(theta[4])*2*(unlist(theta[3])), length(unlist(theta[3])))
  Jac_44 <- as.numeric(unlist(theta[3])^2)
  Jac_4 <- matrix(rbind(c(Jac_41,Jac_42), cbind(Jac_43, Jac_44)), nrow = r34)
  Jac_mat <- matrix(rbind(cbind(Jac_1,Jac_2), cbind(Jac_3, Jac_4)), nrow = r)
  info_reparam1 <- solve(t(Jac_mat)) %*% info %*% solve(Jac_mat)
  rownames(info_reparam1) <- colnames(info_reparam1) <- theta_reparam_names

  Jac_inv_41 <- -as.numeric(unlist(theta[3]))/(2*as.numeric(theta[4]))
  Jac_inv_42 <- diag(1/(2*unlist(theta[3])*as.numeric(theta[4])), length(unlist(theta[3])))
  Jac_inv_43 <- 1
  Jac_inv_44 <- rep(0, length(unlist(theta[3])))
  Jac_inv_4 <- matrix(rbind(cbind(Jac_inv_41, Jac_inv_42), c(Jac_inv_43, Jac_inv_44)), nrow = r34)
  Jac_inv_mat <- matrix(rbind(cbind(Jac_inv_1, Jac_inv_2), cbind(Jac_inv_3, Jac_inv_4)), nrow = r)
  info_reparam2 <- t(Jac_inv_mat) %*% info %*% Jac_inv_mat
  rownames(info_reparam2) <- colnames(info_reparam2) <- theta_reparam_names

  expect_equal(solve(Jac_mat), Jac_inv_mat)
  expect_equal(info_sep, info_reparam1)
  expect_equal(info_reparam1, info_reparam2)

})

