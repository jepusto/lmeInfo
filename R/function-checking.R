
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

check_block_dims <- function(x, m, ni) {
  x_dims <- sapply(x, dim)
  testthat::expect_identical(length(x), m)
  testthat::expect_identical(as.integer(x_dims[1,]), ni)
  testthat::expect_identical(as.integer(x_dims[2,]), ni)
}

test_deriv_dims <- function(mod) {

  vc_est <- extract_varcomp(mod)
  m <- mod$dims$ngrps[[1]]
  ni <- as.integer(table(rev(mod$groups)[[1]]))

  d_Tau <- dV_dreStruct(mod)
  testthat::expect_identical(lengths(d_Tau), lengths(vc_est$Tau))
  # also check that dimensions of d_Tau are consistent with structure of mod

  d_cor <- dV_dcorStruct(mod)
  testthat::expect_identical(length(d_cor), length(vc_est$cor_params))
  # also check that dimensions of d_cor are consistent with structure of mod

  d_var <- dV_dvarStruct(mod)
  testthat::expect_identical(length(d_var), length(vc_est$var_params))
  lapply(d_var, check_block_dims, m = m, ni = ni)
  # also check that dimensions of d_var are consistent with structure of mod

  d_sigma <- build_var_cor_mats(mod, sigma_scale = FALSE)
  # also check that dimensions of d_sigma are consistent with structure of mod

  info_E <- Fisher_info(mod, type = "expected")
  testthat::expect_identical(dim(info_E), rep(length(unlist(vc_est)), 2))

}
