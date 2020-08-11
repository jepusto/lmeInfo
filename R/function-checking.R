
#-------------------------------------------------------------
# Checks that (XWX)^-1 is equivalent to vcov(mod)
#-------------------------------------------------------------

test_Sigma_mats <- function(mod, grps = mod$groups[[1]], sigma_scale = FALSE) {

  Sigma_list <- build_Sigma_mats(mod, invert = TRUE, sigma_scale = sigma_scale)

  # check dimensions
  grp_size <- as.numeric(table(grps))
  dims <- sapply(Sigma_list, dim)
  testthat::expect_equal(grp_size, dims[1,], check.attributes = FALSE)
  testthat::expect_equal(grp_size, dims[2,], check.attributes = FALSE)

  # check that (XWX)^-1 is equivalent to vcov(mod)
  X_design <- model.matrix(mod, data = nlme::getData(mod))
  XVinv <- prod_matrixblock(A = t(X_design), B = Sigma_list, block = attr(Sigma_list, "groups"))
  XWX <- XVinv %*% X_design
  B <- chol2inv(chol(XWX))
  if (!sigma_scale) B <- mod$sigma^2 * B
  testthat::expect_equal(B, vcov(mod), check.attributes = FALSE)
}


#--------------------------------------------------------------------------------
# Checks that dimensions of derivatives are consistent with number of parameters
#--------------------------------------------------------------------------------

expect_correct_block_dims <- function(x, m, ni, is_list = TRUE) {

  if (is_list) {
    return(lapply(x, expect_correct_block_dims, m = m, ni = ni, is_list = FALSE))
  }

  correct_m <- identical(length(x), m)

  x_is_mat <- sapply(x, is.matrix)

  if (all(x_is_mat)) {
    x_dims <- sapply(x, dim)
    correct_dim <- c(
      identical(as.integer(x_dims[1,]), as.integer(ni[names(x)])),
      identical(as.integer(x_dims[2,]), as.integer(ni[names(x)]))
    )
  } else if (all(!x_is_mat)) {
    x_dims <- lengths(x)
    correct_dim <- identical(as.integer(x_dims), as.integer(ni[names(x)]))
  } else {
    correct_dim <- FALSE
  }

  testthat::expect(all(c(correct_m, correct_dim)), "Block dimensions are not correct.")
}

test_deriv_dims <- function(mod) UseMethod("test_deriv_dims")

test_deriv_dims.default <- function(mod) stop("Shouldn't get here!")

test_deriv_dims.gls <- function(mod) {

  vc_est <- extract_varcomp(mod)

  groups <- get_cor_grouping(mod)
  m <- nlevels(groups)
  ni <- table(groups)

  if (!is.null(mod$modelStruct$corStruct)) {
    d_cor <- dV_dcorStruct(mod)
    testthat::expect_identical(length(d_cor), length(vc_est$cor_params))
    expect_correct_block_dims(d_cor, m = m, ni = ni)
  }

  if (!is.null(mod$modelStruct$varStruct)) {
    d_var <- dV_dvarStruct(mod)
    testthat::expect_identical(length(d_var), length(vc_est$var_params))
    expect_correct_block_dims(d_var, m = m, ni = ni)
  }

  d_sigma <- build_var_cor_mats(mod, sigma_scale = FALSE)
  expect_correct_block_dims(d_sigma, m = m, ni = ni, is_list = FALSE)

  info_E <- Fisher_info(mod, type = "expected")
  info_A <- Fisher_info(mod, type = "average")
  r_dim <- rep(length(unlist(vc_est)), 2)

  testthat::expect_identical(dim(info_E), r_dim)
  testthat::expect_identical(dim(info_A), r_dim)
}

test_deriv_dims.lme <- function(mod) {

  vc_est <- extract_varcomp(mod)
  m <- mod$dims$ngrps[names(vc_est$Tau)]
  ni <- lapply(mod$groups[names(vc_est$Tau)], table)
  G <- length(vc_est$Tau)

  if (!is.null(mod$modelStruct$reStruct)) {
    d_Tau <- dV_dreStruct(mod)
    testthat::expect_identical(lengths(d_Tau), lengths(vc_est$Tau))
    mapply(expect_correct_block_dims, x = d_Tau, m = m, ni = ni)
  }

  d_cor <- dV_dcorStruct(mod)
  testthat::expect_identical(length(d_cor), length(vc_est$cor_params))
  expect_correct_block_dims(d_cor, m = m[[G]], ni = ni[[G]])

  d_var <- dV_dvarStruct(mod)
  testthat::expect_identical(length(d_var), length(vc_est$var_params))
  expect_correct_block_dims(d_var, m = m[[G]], ni = ni[[G]])

  d_sigma <- build_var_cor_mats(mod, sigma_scale = FALSE)
  expect_correct_block_dims(d_sigma, m = m[[G]], ni = ni[[G]], is_list = FALSE)

  info_E <- Fisher_info(mod, type = "expected")
  info_A <- Fisher_info(mod, type = "average")
  r_dim <- rep(length(unlist(vc_est)), 2)

  testthat::expect_identical(dim(info_E), r_dim)
  testthat::expect_identical(dim(info_A), r_dim)
}

test_with_FIML <- function(mod) {

  r_dim <- rep(length(unlist(extract_varcomp(mod))), 2)

  dat <- nlme::getData(mod)
  mod_FIML <- suppressWarnings(stats::update(mod, data = dat, method = "ML"))
  mod_FIML$data <- dat

  info_E <- Fisher_info(mod_FIML, type = "expected")
  info_A <- Fisher_info(mod_FIML, type = "average")

  testthat::expect_identical(dim(info_E), r_dim)
  testthat::expect_identical(dim(info_A), r_dim)

}

compare_omit_exclude_complete <- function(mod, dat, NA_vals) {

  dat <- dat
  dat_complete <- dat[!NA_vals,]

  mod_omit <- suppressWarnings(stats::update(mod, data = dat, na.action = "na.omit"))
  mod_exclude <- suppressWarnings(stats::update(mod, data = dat, na.action = "na.exclude"))
  mod_comp <- suppressWarnings(stats::update(mod, data = dat_complete))

  varcomp_omit <- extract_varcomp(mod_omit)
  varcomp_exclude <- extract_varcomp(mod_exclude)
  varcomp_comp <- extract_varcomp(mod_comp)
  testthat::expect_identical(varcomp_omit, varcomp_comp)
  testthat::expect_identical(varcomp_exclude, varcomp_comp)

  EI_omit <- Fisher_info(mod_omit, type = "expected")
  EI_exclude <- Fisher_info(mod_exclude, type = "expected")
  EI_comp <- Fisher_info(mod_comp, type = "expected")
  AI_omit <- Fisher_info(mod_omit, type = "average")
  AI_exclude <- Fisher_info(mod_exclude, type = "average")
  AI_comp <- Fisher_info(mod_comp, type = "average")

  testthat::expect_identical(EI_omit, EI_comp)
  testthat::expect_identical(EI_exclude, EI_comp)
  testthat::expect_identical(AI_omit, AI_comp)
  testthat::expect_identical(AI_exclude, AI_comp)

}

test_after_deleting <- function(mod, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # NA values in response

  dat <- nlme::getData(mod)
  # y_var <- attr(getResponse(mod), "label")
  y_var <- as.character(formula(mod)[[2]])

  NA_vals <- as.logical(rbinom(nrow(dat), size = 1, prob = 0.2))
  dat[[y_var]][NA_vals] <- NA

  compare_omit_exclude_complete(mod, dat, NA_vals)

  # NA values in predictors too

  x_vars <- attr(terms(formula(mod)), "term.labels")
  NA_all <- NA_vals

  for (x in x_vars) {
    NA_vals <- as.logical(rbinom(nrow(dat), size = 1, prob = 0.2 / length(x_vars)))
    dat[[x]][NA_vals] <- NA
    NA_all <- NA_all | NA_vals
  }

  compare_omit_exclude_complete(mod, dat, NA_all)

}

test_after_shuffling <- function(mod, by_var = NULL,
                                 tol_param = 10^-3, tol_info = 10^-3,
                                 test = "info", seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  dat <- nlme::getData(mod)

  if (is.null(by_var)) {
    shuffle <- sample(nrow(dat))
  } else {
    shuffle <- unsplit(tapply(1:nrow(dat), by_var, sample), by_var)
  }
  unshuffle <- order(shuffle)
  dat_shuffle <- dat[shuffle,]

  mod_shuffle <- suppressWarnings(stats::update(mod, data = dat_shuffle))
  mod_shuffle$data <- dat_shuffle

  varcomp_orig <- extract_varcomp(mod)
  varcomp_shuf <- extract_varcomp(mod_shuffle)
  testthat::expect_equal(varcomp_orig, varcomp_shuf, tolerance = tol_param)

  p <- length(unlist(varcomp_orig))
  One <- matrix(1, p, p)
  expected_info_ratio <- Fisher_info(mod, type = "expected") / Fisher_info(mod_shuffle, type = "expected")
  averaged_info_ratio <- Fisher_info(mod, type = "average") / Fisher_info(mod_shuffle, type = "average")

  if (test == "diag-info") {
    testthat::expect_equal(diag(expected_info_ratio), diag(One), tolerance = tol_info, check.attributes = FALSE)
    testthat::expect_equal(diag(averaged_info_ratio), diag(One), tolerance = tol_info, check.attributes = FALSE)
  }

  if (test %in% c("info","full")) {
    testthat::expect_equal(expected_info_ratio, One, tolerance = tol_info, check.attributes = FALSE)
    testthat::expect_equal(averaged_info_ratio, One, tolerance = tol_info, check.attributes = FALSE)
  }

  if (test == "full") {

    unscramble_block <- function(A, unshuffle) {
      A_full <- unblock(A)[unshuffle, unshuffle]
      groups <- attr(A, "groups")[unshuffle]
      A_list <- matrix_list(A_full, fac = groups, dim = "both")
      names(A_list) <- names(A)
      A_list
    }

    R_mat <- build_corr_mats(mod)

    if (!is.null(R_mat)) {
      R_shuff <- unscramble_block(build_corr_mats(mod_shuffle), unshuffle)
      testthat::expect_equal(R_mat, R_shuff, check.attributes = FALSE)
    }

    V_list <- build_var_cor_mats(mod)
    V_shuff <- unscramble_block(build_var_cor_mats(mod_shuffle), unshuffle)
    testthat::expect_equal(V_list, V_shuff, check.attributes = FALSE)

    RE_list <- build_RE_mats(mod)
    RE_shuff <- build_RE_mats(mod_shuffle)
    names(RE_shuff) <- levels(attr(RE_shuff, "groups"))
    RE_shuff <- unscramble_block(RE_shuff, unshuffle)
    testthat::expect_equal(RE_list, RE_shuff, check.attributes = FALSE)

    Sigma_list <- build_Sigma_mats(mod)
    Sigma_shuff <- build_Sigma_mats(mod_shuffle)
    Sigma_shuff <- unscramble_block(Sigma_shuff, unshuffle)
    testthat::expect_equal(Sigma_list, Sigma_shuff, check.attributes = FALSE)

    Tau_params <- unlist(dV_dreStruct(mod), recursive = FALSE)
    cor_params <- dV_dcorStruct(mod)
    if (is.null(cor_params)) cor_params <- list()
    var_params <- dV_dvarStruct(mod)
    if (is.null(var_params)) var_params <- list()
    sigma_sq <- dV_dsigmasq(mod)

    Tau_shuff <- unlist(dV_dreStruct(mod_shuffle), recursive = FALSE)
    Tau_shuff <- lapply(Tau_shuff, unscramble_block, unshuffle = unshuffle)
    cor_shuff <- lapply(dV_dcorStruct(mod_shuffle), unscramble_block, unshuffle = unshuffle)
    var_shuff <- lapply(dV_dvarStruct(mod_shuffle), unscramble_block, unshuffle = unshuffle)
    sigma_sq_shuff <- unscramble_block(dV_dsigmasq(mod_shuffle)[[1]], unshuffle)

    testthat::expect_equal(Tau_params, Tau_shuff, check.attributes = FALSE)
    testthat::expect_equal(cor_params, cor_shuff, check.attributes = FALSE)
    testthat::expect_equal(var_params, var_shuff, check.attributes = FALSE)
    testthat::expect_equal(sigma_sq[[1]], sigma_sq_shuff, check.attributes = FALSE)
  }

}

check_name_order <- function(x_list, group_levels = NULL) {
  if (is.null(group_levels)) group_levels <- levels(attr(x_list, "groups"))
  testthat::expect_identical(names(x_list), group_levels)
}

#--------------------------------------------------------------------
# Checks using REML2
#--------------------------------------------------------------------

check_REML2 <- function(mod, print = FALSE) {
  mod$data <- nlme::getData(mod)
  mod2 <- mod
  mod2$method <- "REML2"

  I_E1 <- Fisher_info(mod, type = "expected")
  I_E2 <- Fisher_info(mod2, type = "expected")
  I_A1 <- Fisher_info(mod, type = "average")
  I_A2 <- Fisher_info(mod2, type = "average")

  if (print) {
    rownames(I_E1) <- rownames(I_E2) <- colnames(I_E1) <- colnames(I_E2) <- NULL
    rownames(I_A1) <- rownames(I_A2) <- colnames(I_A1) <- colnames(I_A2) <- NULL
    print(list(REML_ex = I_E1, REML2_ex = I_E2, REML_av = I_A1, REML2_av = I_A2))
  }

  testthat::expect_equal(I_E1, I_E2)
  testthat::expect_equal(I_A1, I_A2)

}

#--------------------------------------------------------------------
# Checks that lmeInfo::g_mlm() matches scdhlm::g_REML()
#--------------------------------------------------------------------

check_against_scdhlm <- function(mod, p_lmeInfo, r_lmeInfo, p_scdhlm = p_lmeInfo, r_scdhlm, infotype = "expected") {

  g_lmeInfo <- g_mlm(mod, p_lmeInfo, r_lmeInfo)

  g_scdhlm <- scdhlm::g_REML(mod, p_scdhlm, r_scdhlm)

  testthat::expect_equal(g_lmeInfo$p_beta, g_scdhlm$p_beta) # numerator of effect size
  testthat::expect_equal(g_lmeInfo$r_beta, g_scdhlm$r_beta) # squared denominator of effect size
  testthat::expect_equal(g_lmeInfo$delta_AB, g_scdhlm$delta_AB) # unadjusted (REML) effect size estimate
  testthat::expect_equal(g_lmeInfo$nu, g_scdhlm$nu) # degrees of freedom
  testthat::expect_equal(g_lmeInfo$kappa, g_scdhlm$kappa) # constant kappa
  testthat::expect_equal(g_lmeInfo$g_AB, g_scdhlm$g_AB) # corrected effect size estimate
  testthat::expect_equal(g_lmeInfo$SE_g_AB^2, g_scdhlm$V_g_AB) # Approximate variance estimate
  testthat::expect_equal(g_lmeInfo$theta$sigma_sq, g_scdhlm$sigma_sq) # Estimated level-1 variance
  testthat::expect_equal(g_lmeInfo$theta$cor_params, g_scdhlm$phi) # Estimated autocorrelation
  testthat::expect_equal(g_lmeInfo$theta$Tau$case, g_scdhlm$Tau) # Vector of level-2 variance components
  testthat::expect_equal(det(g_lmeInfo$info_inv), det(g_scdhlm$I_E_inv)) # Expected information matrix

}

