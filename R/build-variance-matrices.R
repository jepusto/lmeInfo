
# Construct list of block-diagonal correlation matrices

build_corr_mats <- function(mod) {

  if (is.null(mod$modelStruct$corStruct)) {
    return(NULL)
  } else {
    grps <- nlme::getGroups(mod, form = getGroupsFormula(mod$modelStruct$corStruct))
    # R_list <- as.list(rep(1, nlevels(grps)))
    # names(R_list) <- levels(grps)
    # R_sublist <- nlme::corMatrix(mod$modelStruct$corStruct)
    # R_list[names(R_sublist)] <- R_sublist
    R_list <- nlme::corMatrix(mod$modelStruct$corStruct)
    attr(R_list, "groups") <- grps
    return(R_list)
  }
}

# Construct list of block-diagonal lowest-level var-cov matrices

build_var_cor_mats <- function(mod, R_list = build_corr_mats(mod), sigma_scale = FALSE) {

  sigma_sq <- if (sigma_scale) mod$sigma^2 else 1

  if (is.null(R_list)) {

    # if there is no correlation structure,
    # then build block-diagonals with first available grouping variable

    if (is.null(mod$modelStruct$varStruct)) {
      V_list <- split(rep(sigma_sq, length(mod$groups[[1]])),  mod$groups[[1]])
    } else {
      all_groups <- rev(mod$groups)
      wts <- nlme::varWeights(mod$modelStruct$varStruct)[order(do.call(order, all_groups))]
      V_list <- split(sigma_sq / wts^2, mod$groups[[1]])
    }
    attr(V_list, "groups") <- "diagonal"

  } else {

    # if there is a correlation structure,
    # build block-diagonals according to its grouping structure

    if (is.null(mod$modelStruct$varStruct)) {
      V_list <- if (sigma_scale) lapply(R_list, function(x) x * sigma_sq) else R_list
    } else {
      all_groups <- rev(mod$groups)
      sd_vec <- sqrt(sigma_sq) / nlme::varWeights(mod$modelStruct$varStruct)[order(do.call(order, all_groups))]
      sd_list <- split(sd_vec, attr(R_list, "groups"))
      V_list <- Map(function(R, s) tcrossprod(s) * R, R = R_list, s = sd_list)
    }

    attr(V_list, "groups") <- attr(R_list, "groups")
  }

  return(V_list)
}

# Create block-diagonal covariance structure from Z-design and Tau matrices

ZDZt <- function(D, Z_list) {
  lapply(Z_list, function(z) z %*% D %*% t(z))
}

# Construct list of block-diagonal matrices for each random effects grouping structure

build_RE_mats <- function(mod, sigma_scale = FALSE) {

  # Get random effects structure
  all_groups <- rev(mod$groups)

  if (length(all_groups) == 1) {

    D_mat <- as.matrix(mod$modelStruct$reStruct[[1]])
    if (sigma_scale) D_mat <- mod$sigma^2 * D_mat
    Z_mat <- model.matrix(mod$modelStruct$reStruc, getData(mod))
    row.names(Z_mat) <- NULL
    Z_list <- matrix_list(Z_mat, all_groups[[1]], "row")
    ZDZ_list <- ZDZt(D_mat, Z_list)

    attr(ZDZ_list, "groups") <- all_groups[[1]]

  } else {
    if (sigma_scale) {
      D_list <- lapply(mod$modelStruct$reStruct, function(x) mod$sigma^2 * as.matrix(x))
    } else {
      D_list <- lapply(mod$modelStruct$reStruct, as.matrix)
    }
    Z_mat <- model.matrix(mod$modelStruct$reStruc, getData(mod))
    Z_names <- sapply(strsplit(colnames(Z_mat), ".", fixed=TRUE), function(x) x[1])
    row.names(Z_mat) <- NULL
    Z_levels <- lapply(names(all_groups), function(x) Z_mat[,x==Z_names,drop=FALSE])
    Z_levels <- Map(matrix_list, x = Z_levels, fac = all_groups, dim = "row")

    ZDZ_lists <- Map(ZDZt, D = D_list, Z_list = Z_levels)

    for (i in 2:length(all_groups)) {
      ZDZ_lists[[i]] <- add_bdiag(small_mats = ZDZ_lists[[i-1]],
                                  big_mats = ZDZ_lists[[i]],
                                  crosswalk = all_groups[c(i-1,i)])
    }

    ZDZ_list <- ZDZ_lists[[i]]

    attr(ZDZ_list, "groups") <- all_groups[[i]]

  }

  ZDZ_list

}

build_Sigma_mats <- function(mod, invert = FALSE, sigma_scale = FALSE) {

  if (inherits(mod, "nlme")) stop("not implemented for \"nlme\" objects")

  # lowest-level covariance structure
  V_list <- build_var_cor_mats(mod, sigma_scale = sigma_scale)

  # random effects covariance structure
  ZDZ_list <- build_RE_mats(mod, sigma_scale = sigma_scale)

  V_grps <- attr(V_list, "groups")

  if (identical(V_grps, "diagonal")) {

    Sigma_list <- add_diag_bdiag(V_list, ZDZ_list)
    Sigma_grps <- attr(ZDZ_list, "groups")

  } else {

    # Check if lowest-level covariance structure is nested within RE structure
    ZDZ_grps <- attr(ZDZ_list, "groups")
    group_mapping <- tapply(ZDZ_grps, V_grps, function(x) length(unique(x)))
    nested <- all(group_mapping == 1L)
    if (nested) {
      Sigma_list <- add_bdiag(V_list, ZDZ_list, data.frame(V_grps, ZDZ_grps))
      Sigma_grps <- attr(ZDZ_list, "groups")
    } else {
      V_mat <- unblock(V_list, block = V_grps)
      ZDZ_mat <- unblock(ZDZ_list, block = ZDZ_grps)
      Sigma_list <- V_mat + ZDZ_mat
      Sigma_grps <- factor(rep("A", nrow(Sigma_list)))
    }
  }

  if (invert) {
    Sigma_list <- lapply(Sigma_list, function(x) chol2inv(chol(x)))
  }

  attr(Sigma_list, "groups") <- Sigma_grps

  return(Sigma_list)
}

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

#------------------------------------------------------------------------------
# Create AR(1) correlation and inverse correlation matrices
#------------------------------------------------------------------------------

AR1_corr <- function(phi, times) phi^as.matrix(dist(times))

AR1_corr_block <- function(phi, block, times=NULL) {
  if (is.null(times)) times <- lapply(table(block), seq, from=1)
  lapply(times, AR1_corr, phi=phi)
}


AR1_corr_inv <- function(phi, n)
  diag((rep(1,n) + c(0,rep(phi^2,n-2),0))/(1 - phi^2)) +
  rbind(cbind(rep(0,n-1), diag(rep(-phi/(1 - phi^2),n-1))), rep(0,n)) +
  rbind(rep(0,n), cbind(diag(rep(-phi/(1 - phi^2),n-1)),rep(0,n-1)))

#------------------------------------------------------------------------------
# Create lme covariance and inverse covariance matrices
#------------------------------------------------------------------------------

# create covariance matrix from Tau parameters

Tau_matrix <- function(Tau_coef) {
  q <- (sqrt(1 + 8 * length(Tau_coef)) - 1) / 2
  Tau_mat <- matrix(NA,q,q)
  Tau_mat[upper.tri(Tau_mat, diag=TRUE)] <- Tau_coef
  Tau_mat[lower.tri(Tau_mat)] <- t(Tau_mat)[lower.tri(Tau_mat)]
  return(Tau_mat)
}


# create block-diagonal covariance matrix with AR(1) level-1 error

lmeAR1_cov_block <- function(block, Z_design, theta, times = NULL) {
  if(is.matrix(theta$Tau)) {
    Tau_mat <- theta$Tau
  } else {
    Tau_mat <- Tau_matrix(theta$Tau)
  }
  ZTauZ <- by(Z_design, block, function(z) {
    z_mat <- as.matrix(z)
    z_mat %*% Tau_mat %*% t(z_mat)})
  AR1_mat <- AR1_corr_block(phi=theta$phi, block=block, times=times)
  mapply(function(a,b) a + theta$sigma_sq * b, ZTauZ, AR1_mat, SIMPLIFY = FALSE)
}


# create inverse covariance matrix for a single block

# Note: This function assumes that either
# a) the Tau argument is passed as the eigen-decomposition of Tau
# including an indicator vector for which dimensions to use or
# b) Tau is an invertible matrix.
# Option (a) must be used if Tau is of less than full rank.
# what about 1-dimensional Tau with tau_0^2 = 0?

lmeAR1_cov_inv  <- function(Z_design, sigma_sq, phi, Tau, times=NULL) {
  if (is.null(times)) {
    n <- ifelse(is.vector(Z_design), length(Z_design), dim(Z_design)[1])
    A_inv <- AR1_corr_inv(phi, n) / sigma_sq
  } else {
    A_inv <- solve(AR1_corr(phi, times)) / sigma_sq
  }

  if (is.matrix(Tau)) {
    # otherwise assume that Tau is invertible
    Z_mat <- as.matrix(Z_design)
    Tau_inv <- chol2inv(chol(Tau))
  } else {
    if (sum(Tau$use) == 0) return(A_inv) else {
      # use eigen-decomposition of Tau if available
      Z_mat <- as.matrix(Z_design) %*% Tau$vectors[,Tau$use,drop=FALSE]
      Tau_inv <- diag(1/Tau$values[Tau$use], nrow=sum(Tau$use))
    }
  }
  A_inv_Z <- A_inv %*% Z_mat
  A_inv - A_inv_Z %*% chol2inv(chol(Tau_inv + t(Z_mat) %*% A_inv_Z)) %*% t(A_inv_Z)
}


# create block-diagonal inverse covariance matrix

lmeAR1_cov_block_inv <- function(block, Z_design, theta, times=NULL) {
  Tau_eigen <- eigen(Tau_matrix(theta$Tau), symmetric=TRUE)
  Tau_eigen$use <- round(Tau_eigen$values,14) > 0
  if (is.null(times)) {
    by(Z_design, block, lmeAR1_cov_inv, sigma_sq=theta$sigma_sq, phi=theta$phi, Tau=Tau_eigen)
  } else {
    Z_list <- by(Z_design, block, identity)
    mapply(lmeAR1_cov_inv, Z_design = Z_list, times = times,
           MoreArgs = list(sigma_sq=theta$sigma_sq, phi=theta$phi, Tau=Tau_eigen),
           SIMPLIFY = FALSE)
  }
}
