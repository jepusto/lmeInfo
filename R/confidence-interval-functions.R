
CI_SMD_single <- function(delta, kappa, nu, V_delta, cover, bound) {
  J_nu <- 1 - 3 / (4 * nu - 1)
  start_val <- suppressWarnings(J_nu * (delta + c(-1, 1) * stats::qt(1 - (1 - cover) / 2, df = nu) * sqrt(V_delta)) / kappa)
  L <- kappa * stats::nlminb(start = start_val[1],
                      objective = function(ncp) suppressWarnings((stats::qt((1 - cover) / 2, df = nu, ncp=-ncp) + delta / kappa)^2),
                      lower = -bound, upper = bound)$par
  U <- kappa * stats::nlminb(start = start_val[2],
                      objective = function(ncp) suppressWarnings((stats::qt((1 - cover)  / 2, df = nu, ncp=ncp) - delta / kappa)^2),
                      lower = -bound, upper = bound)$par
  if (L == -bound * kappa) L <- start_val[1] * kappa
  if (U == bound * kappa) U <- start_val[2] * kappa
  c(L,U)
}


## symmetric and approximate non-central t confidence interval ####

#' @title Calculates a confidence interval for a standardized mean difference effect size
#'
#' @description Calculates a confidence interval for a \code{g_mlm} object,
#'   using either a central t distribution (for a symmetric interval) or a
#'   non-central t distribution (for an asymmetric interval).
#'
#' @param g an estimated effect size object of class \code{g_mlm}.
#' @param cover confidence level.
#' @param bound numerical tolerance for non-centrality parameter in
#'   \code{\link[stats]{qt}}.
#' @param symmetric If \code{TRUE} (the default), use a symmetric confidence
#'   interval. If \code{FALSE}, use a non-central t approximation to obtain an
#'   asymmetric confidence interval.
#'
#' @export
#'
#' @return A vector of lower and upper confidence bounds.
#'
#' @examples
#'
#' library(nlme)
#' data(Laski, package = "scdhlm")
#' Laski_RML <- lme(fixed = outcome ~ treatment,
#'                  random = ~ 1 | case,
#'                  correlation = corAR1(0, ~ time | case),
#'                  data = Laski)
#' Laski_g <- g_mlm(Laski_RML, p_const = c(0,1), r_const = c(1,0,1))
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

    suppressWarnings(g_AB + c(-1, 1) * stats::qt(1 - (1 - cover) / 2, df = nu) * SE_g_AB)

  } else {

    CI_SMD_single(delta = delta, kappa = kappa, nu = nu,
                  V_delta = V_delta, cover = cover, bound = bound)
  }

}
