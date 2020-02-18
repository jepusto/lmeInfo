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

Q_matrix <- function(mod) {
  V_inv <- build_Sigma_mats(mod, invert = TRUE, sigma_scale = TRUE)
  X <- model.matrix(mod, data = mod$data)
  Vinv_X <- prod_blockmatrix(V_inv, X)
  XVXinv <- chol2inv(chol(t(X) %*% Vinv_X))
  X_Vinv <- t(Vinv_X)
  VinvX_XVXinv_XVinv <- Vinv_X %*% XVXinv %*% X_Vinv
  block_minus_matrix(V_inv, VinvX_XVXinv_XVinv)
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
#' library(nlme)
#' data(Laski)
#' Laski_RML <- lme(fixed = outcome ~ treatment,
#'                  random = ~ 1 | case,
#'                  correlation = corAR1(0, ~ time | case),
#'                  data = Laski)
#' Fisher_info(Laski_RML, type = "expected")
#'
#' @importFrom stats coef
#' @importFrom stats dist
#' @importFrom stats model.matrix
#' @importFrom stats vcov
#'


Fisher_info <- function(mod, type = "expected") {

  fixed_sigma <- attr(mod$modelStruct, "fixedSigma")

  if (fixed_sigma == TRUE) {
    theta <- extract_varcomp(mod)
    theta <- theta[-length(theta)]
    sigma_sq <- NULL                                       # dV_dsigmasq
  } else {
    theta <- extract_varcomp(mod)
    sigma_sq <- list(build_var_cor_mats(mod, sigma_scale = FALSE))
  }

  theta_names <- vapply(strsplit(names(unlist(theta)), split = "[.]"),
                        function(x) paste(unique(x), collapse = "."), character(1L))
  r <- length(unlist(theta))

  # Calculate derivative matrix-lists

  Tau_params <- dV_dreStruct(mod)                           # random effects structure(s)
  cor_params <- dV_dcorStruct(mod)                          # correlation structure
  var_params <- dV_dvarStruct(mod)                          # variance structure

  # Create a list of derivative matrices
  dV_list <- c(unlist(Tau_params, recursive = FALSE), cor_params, var_params, sigma_sq)

  est_method <- mod$method

  if (type == "expected") {

    if (est_method == "ML") {

      V_inv <- build_Sigma_mats(mod, invert = TRUE, sigma_scale = TRUE)

      # create list with Vinv_dV entries

      Vinv_dV <- lapply(dV_list, prod_blockblock, A = V_inv)

      # calculate I_E

      I_E <- matrix(NA, r, r)

      for (i in 1:r)
        for (j in 1:i)
          I_E[i,j] <- I_E[j,i] <- product_trace_blockblock(Vinv_dV[[i]], Vinv_dV[[j]]) / 2

      rownames(I_E) <- colnames(I_E) <- theta_names
      return(I_E)

    } else if (est_method == "REML") {

      Q_mat <- Q_matrix(mod)

      # create list with Q_dV entries

      Q_dV <- lapply(dV_list, prod_matrixblock, A = Q_mat)

      # calculate I_E

      I_E <- matrix(NA, r, r)

      for (i in 1:r)
        for (j in 1:i)
          I_E[i,j] <- I_E[j,i] <- product_trace(Q_dV[[i]], Q_dV[[j]]) / 2

      rownames(I_E) <- colnames(I_E) <- theta_names
      return(I_E)

    }

  } else if (type == "averaged") {

    V_inv <- build_Sigma_mats(mod, invert = TRUE, sigma_scale = TRUE)

    Vinv_unblock <- unblock(V_inv)

    rhat <- residuals(mod, level = 0) # fixed residuals

    Vinv_rhat <- Vinv_unblock %*% rhat # N*1 matrix

    dVr <- sapply(dV_list, prod_blockmatrix, B = Vinv_rhat, simplify = TRUE) # dV*Vinv*rhat (N * r matrix)

    if (est_method == "ML") {

      I_A <- (t(dVr) %*% prod_blockmatrix(V_inv, dVr)) / 2

      rownames(I_A) <- colnames(I_A) <- theta_names

      return(I_A)

    } else if (est_method == "REML") {

      Q_mat <- Q_matrix(mod)

      I_A <- (t(dVr) %*% Q_mat %*% dVr) / 2

      rownames(I_A) <- colnames(I_A) <- theta_names

      return(I_A)

    }
  }
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
#' library(nlme)
#'
#' Laski_RML <- lme(fixed = outcome ~ treatment,
#'                  random = ~ 1 | case,
#'                  correlation = corAR1(0, ~ time | case),
#'                  data = Laski)
#' varcomp_vcov(Laski_RML, type = "expected")
#'

varcomp_vcov <- function(mod, type = "expected") {

  info_mat <- Fisher_info(mod, type = type)

  chol2inv(chol(info_mat))
}
