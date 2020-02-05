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
  # fill in
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

  theta <- extract_varcomp(mod)

  # Calculate derivative matrix-lists

  Tau_params <- dV_dreStruct(mod)                           # random effects structure(s)
  cor_params <- dV_dcorStruct(mod)                          # correlation structure
  var_params <- dV_dvarStruct(mod)                          # variance structure
  sigma_sq <- build_var_cor_mats(mod, sigma_scale = FALSE)  # sigma^2


  est_method <- mod$method

  if (est_method == "FIML") {

  } else if (est_method == "REML") {

  }

  # # Create list with QdV or (V^-1)dV entries
  # QdV <- rep(1L, r)
  #
  # # calculate I_E
  # r <- sum(lengths(theta))
  #
  # I_E <- matrix(NA, r, r)
  # for (i in 1:r)
  #   for (j in 1:i)
  #     I_E[i,j] <- product_trace(QdV[[i]], QdV[[j]]) / 2
  #
  # I_E[upper.tri(I_E)] <- t(I_E)[upper.tri(I_E)]
  #
  # return(I_E)

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
