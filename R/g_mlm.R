#----------------------------------
# g mlm
#----------------------------------

## estimate adjusted mlm effect size (with associated estimates) for multiple baseline design ####

#' @title Calculates adjusted mlm effect size
#'
#' @description Estimates a standardized mean difference effect size from a
#'   fitted multi-level model, using restricted or full maximum likelihood
#'   methods with small-sample correction, as described in Pustejovsky, Hedges,
#'   & Shadish (2014).
#'
#' @param mod Fitted model of class lmeStruct (estimated using
#'   \code{nlme::lme()}) or of class glsStruct (estimated using
#'   \code{nlme::gls()}), from which to estimate the numerator of the effect
#'   size.
#' @param p_const Vector of constants for calculating numerator of effect size.
#'   Must be the same length as fixed effects in \code{mod}.
#' @param r_const Vector of constants for calculating denominator of effect
#'   size. Must be the same length as the number of variance component
#'   parameters in \code{mod_denom}.
#' @param mod_denom Fitted model of class lmeStruct (estimated using
#'   \code{nlme::lme()}) or of class glsStruct (estimated using
#'   \code{nlme::gls()}), from which to estimate the denominator of the effect
#'   size. If not otherwise specified, the same model will be used for the
#'   numerator and the denominator calculations.
#' @param infotype Type of information matrix. One of \code{"expected"} (the
#'   default), \code{"observed"}, or \code{"average"}.
#' @param separate_variances Logical indicating whether to incorporate separate
#'   level-1 variance components in the calculation of the effect size and
#'   standard error for models with a `varIdent()` variance structure. If
#'   \code{TRUE}, make sure the \code{r_const} matches the parameterization of
#'   the variance component as returned by \code{extract_varcomp(mod,
#'   separate_variances = TRUE)}. Default is \code{FALSE}.
#'
#' @export
#'
#' @return A list with the following components \tabular{ll}{ \code{p_beta} \tab
#'   Numerator of effect size \cr \code{r_theta} \tab Squared denominator of
#'   effect size \cr \code{delta_AB} \tab Unadjusted (mlm) effect size estimate
#'   \cr \code{nu} \tab Estimated denominator degrees of freedom \cr \code{J_nu}
#'   \tab Biased correction factor for effect size estimate \cr \code{kappa}
#'   \tab Scaled standard error of numerator \cr \code{g_AB} \tab Corrected
#'   effect size estimate \cr \code{SE_g_AB} \tab Approximate standard error
#'   estimate \cr \code{cnvg_warn} \tab Indicator that model did not converge
#'   \cr \code{theta} \tab Estimated variance component parameters \cr
#'   \code{info_inv} \tab Inversed information matrix \cr }
#'
#' @references Pustejovsky, J. E., Hedges, L. V., & Shadish, W. R. (2014).
#'   Design-comparable effect sizes in multiple baseline designs: A general
#'   modeling framework. \emph{Journal of Educational and Behavioral Statistics,
#'   39}(4), 211-227. \doi{10.3102/1076998614547577}
#'
#'
#'
#' @examples
#'
#' library(nlme)
#' data(Bryant2016, package = "lmeInfo")
#' Bryant2016_RML1 <- lme(fixed = outcome ~ treatment,
#'                        random = ~ 1 | school/case,
#'                        correlation = corAR1(0, ~ session | school/case),
#'                        data = Bryant2016)
#' Bryant2016_g1 <- g_mlm(Bryant2016_RML1, p_const = c(0,1), r_const = c(1,1,0,1),
#'                        infotype = "expected")
#' print(Bryant2016_g1)
#' summary(Bryant2016_g1)
#'
#'
#' Bryant2016_RML2 <- lme(fixed = outcome ~ treatment,
#'                       random = ~ 1 | school/case,
#'                       correlation = corAR1(0, ~ session | school/case),
#'                       weights = varIdent(form = ~ 1 | treatment),
#'                       data = Bryant2016)
#' Bryant_g <- g_mlm(Bryant2016_RML2, p_const = c(0,1), r_const = c(1,1,0,0,1))
#' Bryant_g_baseline <- g_mlm(Bryant2016_RML2,
#'                            p_const = c(0,1),
#'                            r_const = c(1,1,0,1,0),
#'                            separate_variances = TRUE)
#' Bryant_g_treatment <- g_mlm(Bryant2016_RML2,
#'                             p_const = c(0,1),
#'                             r_const = c(1,1,0,0,1),
#'                             separate_variances = TRUE)
#' print(Bryant_g)
#' print(Bryant_g_baseline)
#' print(Bryant_g_treatment)
#'
#'

g_mlm <- function(mod, p_const, mod_denom = mod, r_const = NULL, infotype = "expected", separate_variances = FALSE) {

  # basic model estimates

  if (inherits(mod, "gls")) {
    beta_coef <- coef(mod)
  } else if (inherits(mod, "lme")) {
    beta_coef <- nlme::fixed.effects(mod)
  } else {
    stop("g_mlm() is only available for lme or gls models. Please specify such a model in the 'mod' argument.")
  }

  if (length(beta_coef) != length(p_const)) {
    stop("The p_const vector must have an entry for every fixed effect coefficient in 'mod'.")
  }

  p_beta <- sum(beta_coef * p_const)
  SE_beta <- sqrt(diag(vcov(mod)))

  if (!inherits(mod_denom, c("gls","lme"))) {
    stop("g_mlm() is only available for lme or gls models. Please specify such a model in the 'mod_denom' argument.")
  }

  theta <- extract_varcomp(mod_denom, separate_variances = separate_variances)

  if (is.null(r_const)) {
    warning("The r_const argument was not specified. Defaulting to r_const equal to all 1's. Are you sure this is right?")
    r_const <- rep(1L, length(unlist(theta)))
  }

  r_theta <- sum(unlist(theta) * r_const)                         # r'theta (sum of var comp)
  delta_AB <- p_beta / sqrt(r_theta)                              # delta_AB
  kappa_sq <- sum(tcrossprod(p_const) * vcov(mod)) / r_theta      # kappa^2
  cnvg_warn <- !is.null(attr(mod_denom,"warning"))                # indicator that RML estimation has not converged

  # calculate inverse Fisher information
  info_inv <- varcomp_vcov(mod_denom, type = infotype, separate_variances = separate_variances)
  SE_theta <- sqrt(diag(info_inv))                                # SE of theta
  nu <- 2 * r_theta^2 / sum(tcrossprod(r_const) * info_inv)       # df
  J_nu <- 1 - 3 / (4 * nu - 1)                                    # bias-correction factor
  g_AB <- J_nu * delta_AB                                         # bias-corrected effect size

  nu_trunc <- max(nu, 2.001)
  J_nu_tr <- 1 - 3 / (4 * nu_trunc - 1)
  V_g_AB <- J_nu^2 * (nu_trunc * kappa_sq / (nu_trunc - 2) + g_AB^2 * (nu_trunc / (nu_trunc - 2) - 1 / J_nu_tr^2))
  SE_g_AB <- sqrt(V_g_AB)

  res <- c(list(p_beta = p_beta, r_theta = r_theta, delta_AB = delta_AB, nu = nu, J_nu = J_nu,
                kappa = sqrt(kappa_sq), g_AB = g_AB, SE_g_AB = SE_g_AB, cnvg_warn=cnvg_warn),
           list(beta = beta_coef, SE_beta = SE_beta),
           list(theta = theta, SE_theta = SE_theta),
           list(info_inv = info_inv, p_const = p_const, r_const = r_const, logLik = mod$logLik))

  class(res) <- "g_mlm"

  return(res)
}

#' @export

summary.g_mlm <- function(object, digits = 3, ...) {

  varcomp <- with(object, cbind(est = c(unlist(theta), "total variance" = r_theta),
                                se = c(unlist(SE_theta), r_theta * sqrt(2 / nu))))

  betas <- with(object, cbind(est = c(beta, "treatment effect at a specified time" = p_beta),
                              se = c(SE_beta, kappa * sqrt(r_theta))))

  ES <- with(object, cbind(est = c("unadjusted effect size" = delta_AB, "adjusted effect size" = g_AB,
                                   "degree of freedom" = nu, "constant kappa" = kappa, logLik = logLik),
                           se = c(SE_g_AB / J_nu, SE_g_AB, NA, NA, NA)))

  print(round(rbind(varcomp, betas, ES), digits), na.print = "")

}

#' @export

print.g_mlm <- function(x, digits = 3, ...) {
  ES <- with(x, cbind(est = c("unadjusted effect size" = delta_AB,
                              "adjusted effect size" = g_AB,
                              "degree of freedom" = nu),
                      se = c(SE_g_AB / J_nu, SE_g_AB, NA)))
  print(round(ES, digits), na.print = "")
}

