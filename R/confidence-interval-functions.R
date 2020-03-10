
CI_SMD_single <- function(delta, kappa, nu, V_delta, cover, bound) {
  J_nu <- 1 - 3 / (4 * nu - 1)
  start_val <- suppressWarnings(J_nu * (delta + c(-1, 1) * qt(1 - (1 - cover) / 2, df = nu) * sqrt(V_delta)) / kappa)
  L <- kappa * nlminb(start = start_val[1],
                      objective = function(ncp) suppressWarnings((qt((1 - cover) / 2, df = nu, ncp=-ncp) + delta / kappa)^2),
                      lower = -bound, upper = bound)$par
  U <- kappa * nlminb(start = start_val[2],
                      objective = function(ncp) suppressWarnings((qt((1 - cover)  / 2, df = nu, ncp=ncp) - delta / kappa)^2),
                      lower = -bound, upper = bound)$par
  if (L == -bound * kappa) L <- start_val[1] * kappa
  if (U == bound * kappa) U <- start_val[2] * kappa
  c(L,U)
}

CI_SMD <- function(delta, kappa, nu, cover = 0.95, bound = 35) {
  J_nu <- 1 - 3 / (4 * nu - 1)
  nu <- pmax(nu, 2.001)
  V_delta = kappa^2 * nu / (nu - 2) + delta^2 * (nu / (nu - 2) - 1 / J_nu^2)
  mapply(CI_SMD_single, delta = delta, kappa = kappa, nu = nu, V_delta = V_delta,
         MoreArgs = list(cover = cover, bound = bound))
}

coverage <- function(delta, CI) CI[1,] < delta & CI[2,] > delta


## symmetric and approximate non-central t confidence interval ####

#' @title Calculates confidence interval for BC-SMD effect size estimates
#'
#' @description Calculates a symmetric confidence interval given a \code{g_mlm} object,
#' based on a central t distribution; and calculates an approximate confidence interval
#' given a \code{g_mlm} object, based on a non-central t approximation.
#'
#' @param g an estimated effect size object of class \code{g_mlm}
#' @param cover confidence level
#' @param bound numerical tolerance for non-centrality parameter in \code{\link{qt}}.
#'
#' @export
#'
#' @return A vector of lower and upper confidence bounds.
#'
#' @examples
#' data(Laski)
#' Laski_RML <- lme(fixed = outcome ~ treatment,
#'                  random = ~ 1 | case,
#'                  correlation = corAR1(0, ~ time | case),
#'                  data = Laski)
#' Laski_g <- g_mlm(Laski_RML, p_const = c(0,1),
#'                   r_const = c(1,0,1), returnModel = FALSE)
#' CI_g(Laski_g, symmetric = TRUE)
#' CI_g(Laski_g, symmetric = FALSE)
#'
#'

CI_g <- function(g, cover = 0.95, bound = 35, symmetric = TRUE) UseMethod("CI_g")

#' @export

CI_g.g_mlm <- function(g, cover = 0.95, bound = 35, symmetric = TRUE) {
  delta <- g$delta_AB
  g_AB <- g$g_AB
  SE_g_AB <- g$SE_g_AB
  kappa <- g$kappa
  nu <- g$nu
  J_nu <- 1 - 3 / (4 * nu - 1)
  V_delta <- kappa^2 * nu / (nu - 2) + delta^2 * (nu / (nu - 2) - 1 / J_nu^2)

  if (symmetric) {

    suppressWarnings(g_AB + c(-1, 1) * qt(1 - (1 - cover) / 2, df = nu) * SE_g_AB)

  } else {

    CI_SMD_single(delta = delta, kappa = kappa, nu = nu,
                  V_delta = V_delta, cover = cover, bound = bound)
  }

}
