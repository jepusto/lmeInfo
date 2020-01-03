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
# Expected Information Matrix
#------------------------------------------------------------------------------

#' @title Calculate expected information matrix
#'
#' @description Calculates the expected information matrix from a fitted linear mixed effects
#' model (lmeStruct object)
#'
#' @param mod Fitted model of class lmeStruct
#'
#' @export
#'
#' @return Expected Information matrix corresponding to variance components of \code{mod}.
#'
#' @examples
#' data(Laski)
#' Laski_RML <- lme(fixed = outcome ~ treatment,
#'                  random = ~ 1 | case,
#'                  correlation = corAR1(0, ~ time | case),
#'                  data = Laski)
#' Info_Expected(Laski_RML)

Info_Expected <- function(mod) {

  theta <- extract_varcomp(mod)
  blocks <- mod$groups
  X_design <- model.matrix(mod, data = mod$data)
  Z_design <- model.matrix(mod$modelStruct$reStruct, data = mod$data)
  times <- attr(mod$modelStruct$corStruct, "covariate")

  Q_mat <- Q_matrix(block, X_design, Z_design, theta, times=times)

  Z_list <- sapply(names(blocks),
                   function(x) Z_design[,grep(x, colnames(Z_design)), drop = FALSE],
                   simplify = FALSE, USE.NAMES = TRUE)


  # Calculate derivative matrix-lists
  sigma_sq <-                                            # sigma^2
  cor_params <- dV_dcorStruct(mod$modelStruct$corStruct)    # correlation structure
  var_params <- dV_dvarStruct(mod$modelStruct$varStruct)    # variance structure
  Tau_params <- mapply(dV_dTau_unstruct,                    # Tau
                       block = blocks, Z_design = Z_list,
                       SIMPLIFY = FALSE)


  # create N * N * r array with QdV entries
  r <- length(unlist(theta))
  QdV <- array(NA, dim = c(dim(Q_mat),r))
  QdV[,,1] <- prod_matrixblock(Q_mat, dV_dsigmasq(block=block, times=times, phi=theta$phi), block=block)
  QdV[,,2] <- prod_matrixblock(Q_mat, dV_dphi(block=block, times=times, phi=theta$phi, sigma_sq=theta$sigma_sq), block=block)
  QdV[,,-2:-1] <- unlist(lapply(dV_dTau_unstruct(block, Z_design), prod_matrixblock, A=Q_mat, block=block))

  # calculate I_E
  I_E <- matrix(NA, r, r)
  for (i in 1:r)
    for (j in 1:i)
      I_E[i,j] <- product_trace(QdV[,,i], QdV[,,j]) / 2
  I_E[upper.tri(I_E)] <- t(I_E)[upper.tri(I_E)]

  return(I_E)

}
