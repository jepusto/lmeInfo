#------------------------------------------------------------------------------
# extract variance components in natural parameterization
#------------------------------------------------------------------------------

#' @title Extract estimated variance components
#'
#' @description Extracts the estimated variance components from a fitted linear
#'   mixed effects model (lmeStruct object) or generalized least squares model
#'   (glsStruct object).
#'
#' @param mod Fitted model of class lmeStruct or glsStruct.
#' @param vector Logical indicating whether to return the variance components as
#'   a numeric vector. Default is \code{FALSE}.
#'
#' @export
#'
#' @return If \code{vector = FALSE}, an object of class \code{varcomp}
#'   consisting of a list of estimated variance components. Models that do not
#'   include correlation structure parameters or variance structure parameters
#'   will have empty lists for those components. If \code{vector = TRUE}, a
#'   numeric vector of estimated variance components.
#'
#' @examples
#'
#' library(nlme)
#' data(Bryant2016)
#' Bryant2016_RML <- lme(fixed = outcome ~ treatment,
#'                       random = ~ 1 | school/case,
#'                       correlation = corAR1(0, ~ session | school/case),
#'                       data = Bryant2016)
#' extract_varcomp(Bryant2016_RML)
#' extract_varcomp(Bryant2016_RML, vector = TRUE)
#'

extract_varcomp <- function(mod, vector) UseMethod("extract_varcomp")

#' @export

extract_varcomp.default <- function(mod, vector = FALSE) {
  mod_class <- paste(class(mod), collapse = "-")
  stop(paste0("Variance components not available for models of class ", mod_class, "."))
}

#' @export

extract_varcomp.gls <- function(mod, vector = FALSE) {

  fixed_sigma <- attr(mod$modelStruct, "fixedSigma")
  sigma_sq <- if (fixed_sigma) NULL else mod$sigma^2                # sigma^2
  cor_params <- as.double(coef(mod$modelStruct$corStruct, FALSE))   # correlation structure
  var_params <- as.double(coef(mod$modelStruct$varStruct, FALSE))   # variance structure

  varcomp <- list(cor_params = cor_params, var_params = var_params, sigma_sq = sigma_sq)

  if (vector) {
    return(unlist(varcomp))
  } else {
    class(varcomp) <- "varcomp"
    return(varcomp)
  }

}

#' @export

extract_varcomp.lme <- function(mod, vector = FALSE) {

  sigma_sq <- mod$sigma^2                                           # sigma^2
  # Tau_params <- coef(mod$modelStruct$reStruct, FALSE) * sigma_sq    # unique coefficients in Tau

  # Calculate Tau_params while taking care of the pdDiag
  RE_params <- coef(mod$modelStruct$reStruct, FALSE)
  Tau_params <- RE_params * RE_params^(as.numeric(grepl("sd", attr(RE_params, "names")))) * sigma_sq
  names(Tau_params) <- mapply(gsub, ".sd", ".var", names(Tau_params), USE.NAMES = FALSE)

  cor_params <- as.double(coef(mod$modelStruct$corStruct, FALSE))   # correlation structure
  var_params <- as.double(coef(mod$modelStruct$varStruct, FALSE))   # variance structure

  # split Tau by grouping variables
  group_names <- names(mod$groups)
  Tau_param_list <- sapply(group_names,
                           function(x) Tau_params[grep(x, names(Tau_params))],
                           simplify = FALSE, USE.NAMES = TRUE)

  fixed_sigma <- attr(mod$modelStruct, "fixedSigma")

  sigma_sq <- if (fixed_sigma) NULL else sigma_sq

  varcomp <- list(Tau = Tau_param_list, cor_params = cor_params, var_params = var_params, sigma_sq = sigma_sq)

  if (vector) {
    return(unlist(varcomp))
  } else {
    class(varcomp) <- "varcomp"
    return(varcomp)
  }

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

#' @title Calculate expected, observed, or average Fisher information matrix
#'
#' @description Calculates the expected, observed, or average Fisher information matrix
#'   from a fitted linear mixed effects model (lmeStruct object)
#'   or generalized least squares model (glsStruct object).
#'
#' @param mod Fitted model of class lmeStruct or glsStruct.
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
#' data(Bryant2016)
#' Bryant2016_RML <- lme(fixed = outcome ~ treatment,
#'                       random = ~ 1 | school/case,
#'                       correlation = corAR1(0, ~ session | school/case),
#'                       data = Bryant2016)
#' Fisher_info(Bryant2016_RML, type = "expected")
#' Fisher_info(Bryant2016_RML, type = "average")
#'
#' @importFrom stats coef
#' @importFrom stats dist
#' @importFrom stats model.matrix
#' @importFrom stats vcov
#' @importFrom stats complete.cases
#' @importFrom stats formula
#' @importFrom stats na.action
#' @importFrom stats rbinom
#' @importFrom stats terms
#'

Fisher_info <- function(mod, type = "expected") {

  theta <- extract_varcomp(mod)
  theta_names <- vapply(strsplit(names(unlist(theta)), split = "[.]"),
                        function(x) paste(unique(x), collapse = "."), character(1L))

  r <- length(unlist(theta))

  # Calculate derivative matrix-lists

  dV_list <- build_dV_list(mod)

  # block-diagonal V^-1
  V_inv <- build_Sigma_mats(mod, invert = TRUE, sigma_scale = TRUE)

  # list with V^-1 dV entries
  Vinv_dV <- lapply(dV_list, prod_blockblock, A = V_inv)

  est_method <- mod$method

  # For REML, need X and M matrices
  if (est_method == "REML") {
    X <- model.matrix(mod, data = nlme::getData(mod))

    # check for columns dropped from model
    col_names <- names(if (inherits(mod, "gls")) coef(mod) else nlme::fixef(mod))
    if (ncol(X) != length(col_names)) X <- X[,col_names,drop=FALSE]

    Vinv_X <- prod_blockmatrix(V_inv, X, block = attr(V_inv, "groups"))
    M <- chol2inv(chol(t(X) %*% Vinv_X))
  }

  if (type == "expected") {

    if (est_method == "ML") {

      # calculate information matrix

      info <- matrix(NA, r, r)

      for (i in 1:r)
        for (j in 1:i)
          info[i,j] <- info[j,i] <- product_trace_blockblock(Vinv_dV[[i]], Vinv_dV[[j]]) / 2

    } else if (est_method == "REML") {

      Vinv_X_M <- Vinv_X %*% M

      # create lists with Xt v^-1 dV entries
      Xt_Vinv_dV <- lapply(Vinv_dV, prod_matrixblock, A = t(X))
      Vinv_dV_Vinv_X_M <- lapply(Vinv_dV, prod_blockmatrix, B = Vinv_X_M)
      dBinv_B <- lapply(Xt_Vinv_dV, function(x) x %*% Vinv_X_M)

      # calculate information matrix

      info <- matrix(NA, r, r)

      for (i in 1:r)
        for (j in 1:i) {
          tr_ij <- product_trace_blockblock(Vinv_dV[[i]], Vinv_dV[[j]]) -
            2 * product_trace(Xt_Vinv_dV[[i]], Vinv_dV_Vinv_X_M[[j]]) +
            product_trace(dBinv_B[[i]], dBinv_B[[j]])
          info[i,j] <- info[j,i] <- tr_ij / 2
        }

    } else if (est_method == "REML2") {

      Q_mat <- Q_matrix(mod)

      # create list with Q_dV entries

      Q_dV <- lapply(dV_list, prod_matrixblock, A = Q_mat)

      # calculate information matrix

      info <- matrix(NA, r, r)

      for (i in 1:r)
        for (j in 1:i)
          info[i,j] <- info[j,i] <- product_trace(Q_dV[[i]], Q_dV[[j]]) / 2

    }

  } else if (type == "average") {

    rhat <- as.matrix(stats::residuals(mod, level = 0)) # fixed residuals
    if (inherits(na.action(mod), "exclude")) rhat <- rhat[-as.integer(na.action(mod)),,drop=FALSE]

    Vinv_rhat <- prod_blockmatrix(A = V_inv, B = rhat)

    dVr <- sapply(dV_list, prod_blockmatrix, B = Vinv_rhat, simplify = TRUE) # dV*Vinv*rhat (N * r matrix)

    if (est_method == "ML") {

      info <- (t(dVr) %*% prod_blockmatrix(V_inv, dVr)) / 2

    } else if (est_method == "REML") {

      Xt_Vinv_dV_Vinv_rhat <- t(Vinv_X) %*% dVr

      info <- (t(dVr) %*% prod_blockmatrix(V_inv, dVr)  -
        t(Xt_Vinv_dV_Vinv_rhat) %*% M %*% Xt_Vinv_dV_Vinv_rhat) / 2

    } else if (est_method == "REML2") {

      Q_mat <- Q_matrix(mod)

      info <- (t(dVr) %*% Q_mat %*% dVr) / 2

    }
  }

  rownames(info) <- colnames(info) <- theta_names
  return(info)

}

#------------------------------------------------------------------------------
# Sampling variance-covariance of variance component parameters
#------------------------------------------------------------------------------

#' @title Estimated sampling variance-covariance of variance component
#'   parameters.
#'
#' @description Estimate the sampling variance-covariance of variance component
#'   parameters from a fitted linear mixed effects model (lmeStruct object)
#'   or generalized least squares model (glsStruct object)
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
#' library(nlme)
#' data(Bryant2016)
#' Bryant2016_RML <- lme(fixed = outcome ~ treatment,
#'                       random = ~ 1 | school/case,
#'                       correlation = corAR1(0, ~ session | school/case),
#'                       data = Bryant2016)
#' varcomp_vcov(Bryant2016_RML, type = "expected")
#'

varcomp_vcov <- function(mod, type = "expected") {

  info_mat <- Fisher_info(mod, type = type)

  res <- chol2inv(chol(info_mat))
  dimnames(res) <- dimnames(info_mat)
  res
}
