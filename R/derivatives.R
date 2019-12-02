#------------------------------------------------------------------------------
# Create first derivative matrices
#------------------------------------------------------------------------------

dV_dsigmasq <- function(block, times, phi)
  AR1_corr_block(phi=phi, block=block, times=times)

dV_dphi <- function(block, times, phi, sigma_sq) {
  if (is.null(times)) times <- lapply(table(block), seq, from=1)
  time_diff <- lapply(times, function(x) as.matrix(dist(x)))
  lapply(time_diff, function(x) sigma_sq * x * phi^pmax(0, x - 1))
}

dV_dTau_index <- function(block, Z_design, tau_index)
  by(Z_design, block, function(Z)
    as.matrix(Z)[,tau_index, drop=FALSE] %*% t(Z)[rev(tau_index),,drop=FALSE])

dV_dTau_unstruct <- function(block, Z_design) {
  Tau_q <- dim(Z_design)[2]
  tau_index <- cbind(unlist(sapply(1:Tau_q, function(x) seq(1,x))),
                     unlist(sapply(1:Tau_q, function(x) rep(x,x))))
  apply(tau_index, 1, function(t) dV_dTau_index(block, Z_design, tau_index=unique(t)))
}
