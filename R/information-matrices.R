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
#' @param separate_variances Logical indicating whether to return the separate
#'   level-1 variance components for each stratum if using \code{varIdent}
#'   function to allow for different variances per stratum. Default is
#'   \code{FALSE}.
#' @param vector Logical indicating whether to return the variance components as
#'   a numeric vector. Default is \code{FALSE}.
#'
#' @export
#'
#' @return If \code{vector = FALSE}, an object of class \code{varcomp} consisting of a
#'   list of estimated variance components. Models that do not include
#'   correlation structure parameters or variance structure parameters will have
#'   empty lists for those components. If \code{vector = TRUE}, a numeric vector
#'   of estimated variance components.
#'
#'   If \code{separate_variances = TRUE} and if \code{weights =
#'   varIdent(form = ~ 1 | Stratum)} is specified in the model fitting, separate
#'   level-1 variance estimates will be returned for each stratum. If
#'   \code{separate_variances = TRUE} but if the weighting structure is not
#'   specified with \code{varIdent}, or if \code{separate_variances = FALSE},
#'   then no separate level-1 variance estimates will be returned.
#'
#' @examples
#'
#' library(nlme)
#' data(Bryant2016)
#' Bryant2016_RML <- lme(fixed = outcome ~ treatment,
#'                       random = ~ 1 | school/case,
#'                       correlation = corAR1(0, ~ session | school/case),
#'                       weights = varIdent(form = ~ 1 | treatment),
#'                       data = Bryant2016)
#' extract_varcomp(Bryant2016_RML, separate_variances = FALSE)
#' extract_varcomp(Bryant2016_RML, separate_variances = TRUE)
#' extract_varcomp(Bryant2016_RML, vector = TRUE)
#'

extract_varcomp <- function(mod, separate_variances, vector) UseMethod("extract_varcomp")

#' @export

extract_varcomp.default <- function(mod, separate_variances = FALSE, vector = FALSE) {
  mod_class <- paste(class(mod), collapse = "-")
  stop(paste0("Variance components not available for models of class ", mod_class, "."))
}

#' @export

extract_varcomp.gls <- function(mod, separate_variances = FALSE, vector = FALSE) {

  fixed_sigma <- attr(mod$modelStruct, "fixedSigma")
  sigma_sq <- if (fixed_sigma) NULL else mod$sigma^2                # sigma^2
  cor_params <- as.double(coef(mod$modelStruct$corStruct, FALSE))   # correlation structure
  var_params <- as.double(coef(mod$modelStruct$varStruct, FALSE))   # variance structure

  varcomp <- list(cor_params = cor_params, var_params = var_params, sigma_sq = sigma_sq)

  # get separate variances when relevant
  if (!is.null(mod$call$weights) && inherits(mod$modelStruct$varStruct, "varIdent") && separate_variances) {
    varStruct <- mod$modelStruct$varStruct
    var_formula <- nlme::getGroupsFormula(varStruct)
    dat <- nlme::getData(mod)
    grps <- stats::model.frame(var_formula, data = dat)
    levels <- levels(grps[,1])
    sigma_sq_grps <- sigma_sq * c(1, var_params^2)
    names(sigma_sq_grps) <- levels
    varcomp$var_params <- NULL
    varcomp$sigma_sq <- sigma_sq_grps
  } else if (separate_variances) {
    warning("The separate_variance option is only relevant for models with a `varIdent()` variance structure.")
  }

  if (vector) {
    return(unlist(varcomp))
  } else {
    class(varcomp) <- "varcomp"
    return(varcomp)
  }

}

#' @export

extract_varcomp.lme <- function(mod, separate_variances = FALSE, vector = FALSE) {

  sigma_sq <- mod$sigma^2                                           # sigma^2
  # Tau_params <- coef(mod$modelStruct$reStruct, FALSE) * sigma_sq    # unique coefficients in Tau

  # Calculate Tau_params while taking care of the pdDiag
  RE_params <- coef(mod$modelStruct$reStruct, FALSE)
  Tau_params <- RE_params * RE_params^(as.numeric(grepl("sd", attr(RE_params, "names")))) * sigma_sq
  names(Tau_params) <- mapply(gsub, ".sd", ".var", names(Tau_params), USE.NAMES = FALSE)

  # split Tau by grouping variables
  group_names <- lapply(mod$modelStruct$reStruct, function(x) attr(x, "Dimnames")[[1]])
  group_regx <- paste0(names(group_names), ".+\\(",lapply(group_names, paste, collapse = "|"))
  names(group_regx) <- names(group_names)

  Tau_param_list <- sapply(group_regx,
                           function(x) Tau_params[grepl(x, names(Tau_params))],
                           simplify = FALSE, USE.NAMES = TRUE)

  cor_params <- as.double(coef(mod$modelStruct$corStruct, FALSE))   # correlation structure
  var_params <- as.double(coef(mod$modelStruct$varStruct, FALSE))   # variance structure

  fixed_sigma <- attr(mod$modelStruct, "fixedSigma")

  sigma_sq <- if (fixed_sigma) NULL else sigma_sq

  varcomp <- list(Tau = Tau_param_list, cor_params = cor_params, var_params = var_params, sigma_sq = sigma_sq)

  # get separate variances when relevant
  if (!is.null(mod$call$weights) && inherits(mod$modelStruct$varStruct, "varIdent") && separate_variances) {
    varStruct <- mod$modelStruct$varStruct
    var_formula <- nlme::getGroupsFormula(varStruct)
    dat <- nlme::getData(mod)
    grps <- stats::model.frame(var_formula, data = dat)
    levels <- unique(grps[,1])
    sigma_sq_grps <- sigma_sq * c(1, var_params^2)
    names(sigma_sq_grps) <- levels
    varcomp$var_params <- NULL
    varcomp$sigma_sq <- sigma_sq_grps
  } else if (separate_variances) {
    warning("The separate_variance option is only relevant for models with a `varIdent()` variance structure.")
  }

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
#' @description Calculates the expected, observed, or average Fisher information
#'   matrix from a fitted linear mixed effects model (lmeStruct object) or
#'   generalized least squares model (glsStruct object).
#'
#' @param mod Fitted model of class lmeStruct or glsStruct.
#' @param type Type of information matrix. One of \code{"expected"} (the
#'   default), \code{"observed"}, or \code{"average"}.
#' @param separate_variances Logical indicating whether to return the Fisher
#'   information matrix for separate level-1 variance components if using
#'   \code{varIdent} function to allow for different variances per stratum.
#'   Default is \code{FALSE}.
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
#' Bryant2016_RML2 <- lme(fixed = outcome ~ treatment,
#'                       random = ~ 1 | school/case,
#'                       correlation = corAR1(0, ~ session | school/case),
#'                       weights = varIdent(form = ~ 1 | treatment),
#'                       data = Bryant2016)
#' Fisher_info(Bryant2016_RML2, separate_variances = TRUE)
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

Fisher_info <- function(mod, type = "expected", separate_variances = FALSE) {

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

  # info mat for seperate variances
  if (!is.null(mod$call$weights) & separate_variances) {
    varFunc <- sub("\\(.*", "", mod$call$weights)[1]
    if (varFunc == "varIdent") {
      theta_reparam <- extract_varcomp(mod, separate_variances = TRUE)
      theta_reparam_names <- vapply(strsplit(names(unlist(theta_reparam)), split = "[.]"),
                            function(x) paste(unique(x), collapse = "."), character(1L))
      r12 <- length(unlist(theta[c(1,2)]))
      r34 <- length(unlist(theta[c(3,4)]))
      r <- length(unlist(theta))
      Jac_1 <- diag(1, r12)
      Jac_2 <- matrix(0, nrow = r12, ncol = r34)
      Jac_3 <- t(Jac_2)
      Jac_41 <- rep(0, length(unlist(theta[3])))
      Jac_42 <- 1
      Jac_43 <- diag(as.numeric(theta[4])*2*(unlist(theta[3])), length(unlist(theta[3])))
      Jac_44 <- as.numeric(unlist(theta[3])^2)
      Jac_4 <- matrix(rbind(c(Jac_41,Jac_42), cbind(Jac_43, Jac_44)), nrow = r34)
      Jac_mat <- matrix(rbind(cbind(Jac_1,Jac_2), cbind(Jac_3, Jac_4)), nrow = r)
      info <- solve(t(Jac_mat)) %*% info %*% solve(Jac_mat)
      rownames(info) <- colnames(info) <- theta_reparam_names
    } else {
      warning("The `Fisher_info()` returns information matrix for the variance components that includes the separate level-1 variances only when the variance structure is specified with `varIdent()` in the model.")
    }
  } else if (is.null(mod$call$weights) & separate_variances) {
    warning("The separate_variance option is only relevant when the variance structure is specified with `varIdent()`.")
  }

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
#' Bryant2016_RML2 <- lme(fixed = outcome ~ treatment,
#'                       random = ~ 1 | school/case,
#'                       correlation = corAR1(0, ~ session | school/case),
#'                       weights = varIdent(form = ~ 1 | treatment),
#'                       data = Bryant2016)
#' varcomp_vcov(Bryant2016_RML, separate_variances = TRUE)
#'

varcomp_vcov <- function(mod, type = "expected", separate_variances = FALSE) {

  if (separate_variances) {
    info_mat <- Fisher_info(mod, type = type, separate_variances = TRUE)
  } else {
    info_mat <- Fisher_info(mod, type = type, separate_variances = FALSE)
  }

  res <- chol2inv(chol(info_mat))
  dimnames(res) <- dimnames(info_mat)
  res
}
