#------------------------------------------------------------------------------
# extract variance components in natural parameterization
#------------------------------------------------------------------------------

extract_varcomp <- function(mod) {

  sigma_sq <- mod$sigma^2                                           # sigma^2
  cor_params <- as.double(coef(mod$modelStruct$corStruct, FALSE))   # correlation structure
  var_params <- as.double(coef(mod$modelStruct$varStruct, FALSE))   # variance structure
  Tau_params <- coef(mod$modelStruct$reStruct, FALSE) * sigma_sq    # unique coefficients in Tau

  # split Tau by grouping variables
  group_names <- names(mod$groups)
  Tau_param_list <- sapply(group_names,
                           function(x) Tau_params[grep(x, names(Tau_params))],
                           simplify = FALSE, USE.NAMES = TRUE)

  varcomp <- list(Tau = Tau_param_list, cor_params = cor_params, var_params = var_params, sigma_sq=sigma_sq)

  class(varcomp) <- "varcomp"
  return(varcomp)
}

#------------------------------------------------------------------------------
# create Q matrix
#------------------------------------------------------------------------------

Q_matrix <- function(block, X_design, Z_design, theta, times=NULL) {
  V_inv <- lmeAR1_cov_block_inv(block=block, Z_design=Z_design, theta=theta, times=times)
  Vinv_X <- prod_blockmatrix(V_inv, X_design, block = block)
  VinvX_XVXinv_XVinv <- Vinv_X %*% chol2inv(chol(t(X_design) %*% Vinv_X)) %*% t(Vinv_X)
  block_minus_matrix(V_inv, VinvX_XVXinv_XVinv, block)
}

#------------------------------------------------------------------------------
# Information Matrices
#------------------------------------------------------------------------------

#' @title Calculate expected, observed, or average information matrix
#'
#' @description Calculates the expected, observed, or average information matrix
#'   from a fitted linear mixed effects model (lmeStruct object)
#'
#' @param mod Fitted model of class lmeStruct
#' @param type Type of information matrix. One of \code{"expected"} (the default),
#'   \code{"observed"}, or \code{"average"}.
#'
#' @export
#'
#' @return Information matrix corresponding to variance component parameters of
#'   \code{mod}.
#'
#' @examples
#'
#' data(Laski)
#' Laski_RML <- lme(fixed = outcome ~ treatment,
#'                  random = ~ 1 | case,
#'                  correlation = corAR1(0, ~ time | case),
#'                  data = Laski)
#' Fisher_info(Laski_RML, type = "expected")
#'

Fisher_info <- function(mod, type = "expected") {

  theta <- extract_varcomp(mod)

  # Calculate derivative matrix-lists

  Tau_params <- dV_dreStruct(mod)                           # random effects structure(s)
  cor_params <- dV_dcorStruct(mod)                          # correlation structure
  var_params <- dV_dvarStruct(mod)                          # variance structure
  sigma_sq <- build_var_cor_mats(mod, sigma_scale = FALSE)  # sigma^2


  # Create list with QdV or (V^-1)dV entries

  # calculate I_E

  I_E <- matrix(NA, r, r)
  for (i in 1:r)
    for (j in 1:i)
      I_E[i,j] <- product_trace(QdV[,,i], QdV[,,j]) / 2
  I_E[upper.tri(I_E)] <- t(I_E)[upper.tri(I_E)]

  return(I_E)

}


#------------------------------------------------------------------------------
# Sampling variance-covariance of variance component parameters
#------------------------------------------------------------------------------

#' @title Estimated sampling variance-covariance of variance component
#'   parameters.
#'
#' @description Estimate the sampling variance-covariance of variance component
#'   parameters from a fitted linear mixed effects model (lmeStruct object)
#'   using the inverse Fisher information.
#'
#' @inheritParams Fisher_info
#'
#' @export
#'
#' @return Sampling variance-covariance matrix corresponding to variance
#'   component parameters of \code{mod}.
#'
#' @examples
#'
#' data(Laski)
#' Laski_RML <- lme(fixed = outcome ~ treatment,
#'                  random = ~ 1 | case,
#'                  correlation = corAR1(0, ~ time | case),
#'                  data = Laski)
#' varcomp_vcov(Laski_RML, type = "expected")
#'

varcomp_vcov <- function(mod, type = "expected") {

  Info <- Fisher_info(mod)

  chol2inv(chol(InfO))
}
