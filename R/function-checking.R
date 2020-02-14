
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

  mod_FIML <- update(mod, method = "ML")

  info_E <- Fisher_info(mod, type = "expected")
  info_A <- Fisher_info(mod, type = "averaged")

  testthat::expect_identical(dim(info_E), r_dim)
  testthat::expect_identical(dim(info_A), r_dim)

}
