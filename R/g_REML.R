g_REML <- function(mod, p_const, r_const, infotype = "expected", returnModel = TRUE) {

  # basic model estimates
  p_beta <- sum(nlme::fixed.effects(mod) * p_const)               # p'Beta
  theta <- extract_varcomp(mod)                                   # full theta vector
  r_theta <- sum(unlist(theta) * r_const)                         # r'theta (sum of var comp)
  delta_AB <- p_beta / sqrt(r_theta)                              # delta_AB
  kappa_sq <- (t(p_const) %*% vcov(mod) %*% p_const) / r_theta    # kappa^2
  cnvg_warn <- !is.null(attr(mod,"warning"))                      # indicator that RML estimation has not converged

  # calculate inverse Fisher information
  info <- Fisher_info(mod, type = infotype)
  info_inv <- chol2inv(chol(info))
  rownames(info_inv) <- colnames(info_inv) <- rownames(info)
  SE_theta <- sqrt(diag(info_inv))                                # SE of theta

  nu <- 2 * r_theta^2 / (t(r_const) %*% info_inv %*% r_const)      # df
  J_nu <- 1 - 3 / (4 * nu - 1)                                    # bias-correction factor
  g_AB <- J_nu * delta_AB                                         # bias-corrected effect size

  nu_trunc <- max(nu, 2.001)
  J_nu_tr <- 1 - 3 / (4 * nu_trunc - 1)
  V_g_AB <- J_nu^2 * (nu_trunc * kappa_sq / (nu_trunc - 2) + g_AB^2 * (nu_trunc / (nu_trunc - 2) - 1 / J_nu_tr^2))
  SE_g_AB <- sqrt(V_g_AB)

  res <- c(list(p_beta = p_beta, r_theta = r_theta, delta_AB = delta_AB, nu = nu, J_nu = J_nu,
                kappa = sqrt(kappa_sq), g_AB = g_AB, SE_g_AB = SE_g_AB, cnvg_warn=cnvg_warn),
           list(theta = theta, SE_theta = SE_theta),
           list(info_inv = info_inv, p_const = p_const, r_const = r_const))

  if (returnModel) {
    res <- c(res, mod)
  }

  class(res) <- "g_REML"

  return(res)
}

summary.g_REML <- function(object, ...) {
  varcomp <- with(object, cbind(est = c(unlist(theta), "total variance" = r_theta),
                                se = c(unlist(SE_theta), r_theta * sqrt(2 / nu))))
  betas <- with(object, cbind(est = c(coefficients$fixed, "treatment effect at a specified time" = p_beta),
                              se = c(sqrt(diag(varFix)), kappa * sqrt(r_theta))))
  ES <- with(object, cbind(est = c("unadjusted effect size" = delta_AB, "adjusted effect size" = g_AB,
                                   "degree of freedom" = nu, "constant kappa" = kappa, logLik = logLik),
                           se = c(SE_g_AB / J_nu, SE_g_AB, NA, NA, NA)))

  round(rbind(varcomp, betas, ES), 3)
}

print.g_REML <- function(object, ...) {
  ES <- with(object, cbind(est = c("unadjusted effect size" = delta_AB,
                                   "adjusted effect size" = g_AB,
                                   "degree of freedom" = nu),
                           se = c(SE_g_AB / J_nu, SE_g_AB, NA)))
  round(ES, 3)
}

