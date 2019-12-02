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
