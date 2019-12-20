
# from clubSandwich

ZDZt <- function(D, Z_list) {
  lapply(Z_list, function(z) z %*% D %*% t(z))
}

targetVariance <- function(obj, cluster = nlme::getGroups(obj, level = 1)) {

  if (inherits(obj, "nlme")) stop("not implemented for \"nlme\" objects")

  all_groups <- rev(obj$groups)
  smallest_groups <- all_groups[[1]]
  largest_groups <- all_groups[[length(all_groups)]]

  # Get level-1 variance-covariance structure as V_list

  if (is.null(obj$modelStruct$corStruct)) {
    if (is.null(obj$modelStruct$varStruct)) {
      V_list <- matrix_list(rep(1, length(smallest_groups)), smallest_groups, "both")
    } else {
      wts <- nlme::varWeights(obj$modelStruct$varStruct)[order(do.call(order, all_groups))]
      V_list <- matrix_list(1 / wts^2, smallest_groups, "both")
    }
  } else {
    R_list <- as.list(rep(1, nlevels(smallest_groups)))
    names(R_list) <- levels(smallest_groups)
    R_sublist <- nlme::corMatrix(obj$modelStruct$corStruct)
    R_list[names(R_sublist)] <- R_sublist
    if (is.null(obj$modelStruct$varStruct)) {
      V_list <- R_list
    } else {
      sd_vec <- 1 / nlme::varWeights(obj$modelStruct$varStruct)[order(do.call(order, all_groups))]
      sd_list <- split(sd_vec, smallest_groups)
      V_list <- Map(function(R, s) tcrossprod(s) * R, R = R_list, s = sd_list)
    }
  }

  # Get random effects structure

  if (length(all_groups) == 1) {
    D_mat <- as.matrix(obj$modelStruct$reStruct[[1]])
    Z_mat <- model.matrix(obj$modelStruct$reStruc, getData(obj))
    row.names(Z_mat) <- NULL
    Z_list <- matrix_list(Z_mat, all_groups[[1]], "row")
    ZDZ_list <- ZDZt(D_mat, Z_list)
    target_list <- Map("+", ZDZ_list, V_list)
  } else {
    D_list <- lapply(obj$modelStruct$reStruct, as.matrix)
    Z_mat <- model.matrix(obj$modelStruct$reStruc, getData(obj))
    Z_names <- sapply(strsplit(colnames(Z_mat), ".", fixed=TRUE), function(x) x[1])
    row.names(Z_mat) <- NULL
    Z_levels <- lapply(names(all_groups), function(x) Z_mat[,x==Z_names,drop=FALSE])
    Z_levels <- Map(matrix_list, x = Z_levels, fac = all_groups, dim = "row")
    ZDZ_lists <- Map(ZDZt, D = D_list, Z_list = Z_levels)
    ZDZ_lists[[1]] <- Map("+", ZDZ_lists[[1]], V_list)
    for (i in 2:length(all_groups)) {
      ZDZ_lists[[i]] <- add_bdiag(small_mats = ZDZ_lists[[i-1]],
                                  big_mats = ZDZ_lists[[i]],
                                  crosswalk = all_groups[c(i-1,i)])
    }
    target_list <- ZDZ_lists[[i]]
  }

  # check if clustering level is higher than highest level of random effects

  tb_groups <- table(largest_groups)
  tb_cluster <- table(cluster)

  if (length(tb_groups) < length(tb_cluster) |
      any(as.vector(tb_groups) != rep(as.vector(tb_cluster), length.out = length(tb_groups))) |
      any(names(tb_groups) != rep(names(tb_cluster), length.out = length(tb_groups)))) {

    # check that random effects are nested within clusters
    tb_cross <- table(largest_groups, cluster)
    nested <- apply(tb_cross, 1, function(x) sum(x > 0) == 1)
    if (!all(nested)) stop("Random effects are not nested within clustering variable.")

    # expand target_list to level of clustering
    crosswalk <- data.frame(largest_groups, cluster)
    target_list <- add_bdiag(small_mats = target_list,
                             big_mats = matrix_list(rep(0, length(cluster)), cluster, dim = "both"),
                             crosswalk = crosswalk)
  }

  return(target_list)
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
