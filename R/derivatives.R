#------------------------------------------------------------------------------
# Create first derivative matrices
#------------------------------------------------------------------------------



# Derivative of unstructured random effects covariance matrix

dV_dTau_index <- function(tau_index, Z_blocks) {
  lapply(Z_blocks, function(Z)
    Z[,tau_index, drop=FALSE] %*% t(Z)[rev(tau_index),,drop=FALSE])
}

dV_dTau_unstruct <- function(block, Z_design) {
  Tau_q <- dim(Z_design)[2]
  Z_blocks <- by(Z_design, block, as.matrix)
  tau_index <- cbind(unlist(sapply(1:Tau_q, function(x) seq(1,x))),
                     unlist(sapply(1:Tau_q, function(x) rep(x,x))))
  apply(tau_index, 1, function(t) dV_dTau_index(unique(t), Z_blocks = Z_blocks))
}
