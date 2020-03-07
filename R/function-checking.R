
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
  X_design <- model.matrix(mod, data = mod$data)
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
      identical(as.integer(x_dims[1,]), as.integer(ni)),
      identical(as.integer(x_dims[2,]), as.integer(ni))
    )
  } else if (all(!x_is_mat)) {
    x_dims <- lengths(x)
    correct_dim <- identical(as.integer(x_dims), as.integer(ni))
  } else {
    correct_dim <- FALSE
  }

  testthat::expect(all(c(correct_m, correct_dim)), "Block dimensions are not correct.")
}

test_deriv_dims <- function(mod) {

  vc_est <- extract_varcomp(mod)
  m <- mod$dims$ngrps[names(vc_est$Tau)]
  ni <- lapply(mod$groups[names(vc_est$Tau)], table)
  G <- length(vc_est$Tau)

  d_Tau <- dV_dreStruct(mod)
  testthat::expect_identical(lengths(d_Tau), lengths(vc_est$Tau))
  mapply(expect_correct_block_dims, x = d_Tau, m = m, ni = ni)

  d_cor <- dV_dcorStruct(mod)
  testthat::expect_identical(length(d_cor), length(vc_est$cor_params))
  expect_correct_block_dims(d_cor, m = m[[G]], ni = ni[[G]])

  d_var <- dV_dvarStruct(mod)
  testthat::expect_identical(length(d_var), length(vc_est$var_params))
  expect_correct_block_dims(d_var, m = m[[G]], ni = ni[[G]])

  d_sigma <- build_var_cor_mats(mod, sigma_scale = FALSE)
  expect_correct_block_dims(d_sigma, m = m[[G]], ni = ni[[G]], is_list = FALSE)

  info_E <- Fisher_info(mod, type = "expected")
  info_A <- Fisher_info(mod, type = "averaged")
  r_dim <- rep(length(unlist(vc_est)), 2)

  testthat::expect_identical(dim(info_E), r_dim)
  testthat::expect_identical(dim(info_A), r_dim)
}

test_with_FIML <- function(mod) {

  r_dim <- rep(length(unlist(extract_varcomp(mod))), 2)

  mod_FIML <- suppressWarnings(stats::update(mod, data = nlme::getData(mod), method = "ML"))

  info_E <- Fisher_info(mod_FIML, type = "expected")
  info_A <- Fisher_info(mod_FIML, type = "averaged")

  testthat::expect_identical(dim(info_E), r_dim)
  testthat::expect_identical(dim(info_A), r_dim)

}

test_after_shuffling <- function(mod, tol_param = 10^-5, tol_info = .03) {

  dat <- nlme::getData(mod)
  shuffle <- sample(nrow(dat))
  dat_shuffle <- dat[shuffle,]

  mod_shuffle <- stats::update(mod, data = dat_shuffle)

  varcomp_orig <- extract_varcomp(mod)
  varcomp_shuf <- extract_varcomp(mod_shuffle)
  testthat::expect_equal(varcomp_orig, varcomp_shuf, tolerance = tol_param)

  p <- length(unlist(varcomp_orig))
  One <- matrix(1, p, p)
  expected_info_ratio <- Fisher_info(mod, type = "expected") / Fisher_info(mod_shuffle, type = "expected")
  averaged_info_ratio <- Fisher_info(mod, type = "averaged") / Fisher_info(mod_shuffle, type = "averaged")
  testthat::expect_equal(expected_info_ratio, One, tolerance = tol_info, check.attributes = FALSE)
  testthat::expect_equal(averaged_info_ratio, One, tolerance = tol_info, check.attributes = FALSE)

  # extract_varcomp(mod) %>% unlist()
  # extract_varcomp(mod_shuffle) %>% unlist()
  #
  # R_mat <- unblock(build_corr_mats(mod))[shuffle, shuffle]
  # R_shuff <- unblock(build_corr_mats(mod_shuffle))
  # expect_equal(R_mat, R_shuff)

}

check_name_order <- function(x_list, group_levels = NULL) {
  if (is.null(group_levels)) group_levels <- levels(attr(x_list, "groups"))
  testthat::expect_identical(names(x_list), group_levels)
}

build_block_matrices <- function(mod) {

  R_list <- build_corr_mats(mod)

  check_name_order(R_list)

  all_groups <- rev(mod$groups)

  if (!is.null(mod$modelStruct$varStruct)) {
    sd_vec <- mod$sigma / nlme::varWeights(mod$modelStruct$varStruct)[order(do.call(order, all_groups))]
    sd_list <- split(sd_vec, attr(R_list, "groups"))

    check_name_order(sd_list, group_levels = levels(attr(R_list, "groups")))
  }

  V_list <- build_var_cor_mats(mod, sigma_scale = TRUE)

  check_name_order(V_list)

  ZDZ_list <- build_RE_mats(mod, sigma_scale = TRUE)
  check_name_order(ZDZ_list)

  Tau_params <- dV_dreStruct(mod)
  sapply(Tau_params, check_name_order)

  cor_params <- dV_dcorStruct(mod)
  sapply(cor_params, check_name_order)

  var_params <- dV_dvarStruct(mod)
  sapply(var_params, check_name_order)

}

#--------------------------------------------------------------------
# Checks using REML2
#--------------------------------------------------------------------

check_REML2 <- function(mod, print = FALSE) {
  mod2 <- mod
  mod2$method <- "REML2"

  I_E1 <- Fisher_info(mod, type = "expected")
  I_E2 <- Fisher_info(mod2, type = "expected")
  I_A1 <- Fisher_info(mod, type = "averaged")
  I_A2 <- Fisher_info(mod2, type = "averaged")

  if (print) {
    rownames(I_E1) <- rownames(I_E2) <- colnames(I_E1) <- colnames(I_E2) <- NULL
    rownames(I_A1) <- rownames(I_A2) <- colnames(I_A1) <- colnames(I_A2) <- NULL
    print(list(REML_ex = I_E1, REML2_ex = I_E2, REML_av = I_A1, REML2_av = I_A2))
  }

  testthat::expect_equal(I_E1, I_E2)
  testthat::expect_equal(I_A1, I_A2)

}

#--------------------------------------------------------------------
# Checks that lmeInfo::g_REML() matches scdhlm::g_REML()
#--------------------------------------------------------------------

check_against_scdhlm <- function(mod, p_lmeInfo, r_lmeInfo, p_scdhlm = p_lmeInfo, r_scdhlm, infotype = "expected") {

  g_lmeInfo <- lmeInfo::g_REML(mod, p_lmeInfo, r_lmeInfo)

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
