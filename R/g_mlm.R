#----------------------------------
# g mlm
#----------------------------------

## estimate adjusted mlm effect size (with associated estimates) for multiple baseline design ####

#' @title Calculates adjusted mlm effect size
#'
#' @description Estimates a design-comparable standardized mean difference effect size based on data
#' from a multiple baseline design, using adjusted mlm method as described in Pustejovsky, Hedges,
#' & Shadish (2014). Note that the data must contain one row per measurement occasion per case.
#'
#' @param mod Fitted model of class lmeStruct (estimated using \code{nlme::lme()})
#' @param p_const Vector of constants for calculating numerator of effect size.
#' Must be the same length as fixed effects in \code{mod}.
#' @param r_const Vector of constants for calculating denominator of effect size.
#' Must be the same length as the number of variance component parameters in \code{mod}.
#' @param infotype Type of information matrix. One of \code{"expected"} (the default),
#'   \code{"observed"}, or \code{"average"}.
#' @param returnModel (Optional) If true, the fitted input model is included in the return.
#'
#' @export
#'
#' @return A list with the following components
#' \tabular{ll}{
#' \code{p_beta} \tab Numerator of effect size \cr
#' \code{r_theta} \tab Squared denominator of effect size \cr
#' \code{delta_AB} \tab Unadjusted (mlm) effect size estimate \cr
#' \code{nu} \tab Estimated denominator degrees of freedom \cr
#' \code{J_nu} \tab Biased correction factor for effect size estimate \cr
#' \code{kappa} \tab Scaled standard error of numerator \cr
#' \code{g_AB} \tab Corrected effect size estimate \cr
#' \code{SE_g_AB} \tab Approximate standard error estimate \cr
#' \code{cnvg_warn} \tab Indicator that model did not converge \cr
#' \code{theta} \tab Estimated variance component parameters \cr
#' \code{info_inv} \tab Inversed information matrix \cr
#' }
#'
#' @references Pustejovsky, J. E., Hedges, L. V., & Shadish, W. R. (2014).
#' Design-comparable effect sizes in multiple baseline designs: A general modeling framework.
#' \emph{Journal of Educational and Behavioral Statistics, 39}(4), 211-227. doi:\href{http://doi.org/10.3102/1076998614547577}{10.3102/1076998614547577}
#'
#' @examples
#'
#' library(nlme)
#' data(Laski, package = "scdhlm")
#' Laski_RML <- lme(fixed = outcome ~ treatment,
#'                  random = ~ 1 | case,
#'                  correlation = corAR1(0, ~ time | case),
#'                  data = Laski)
#' summary(Laski_RML)
#' g_mlm(Laski_RML, p_const = c(0,1), r_const = c(1,0,1), returnModel = FALSE)
#'
#' data(Schutte, package = "scdhlm")
#' Schutte$trt.week <- with(Schutte, unlist(tapply((treatment=="treatment") * week,
#'          list(treatment,case), function(x) x - min(x))) + (treatment=="treatment"))
#' Schutte$week <- Schutte$week - 9
#' Schutte_RML <- lme(fixed = fatigue ~ week + treatment + trt.week,
#'                    random = ~ week | case,
#'                    correlation = corAR1(0.2, ~ week | case),
#'                    data = subset(Schutte, case != 4))
#' summary(Schutte_RML)
#' Schutte_g <- g_mlm(Schutte_RML, p_const = c(0,0,1,7), r_const = c(1,0,0,0,1))
#' summary(Schutte_g)

g_mlm <- function(mod, p_const, r_const, infotype = "expected", returnModel = TRUE) {

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

  class(res) <- "g_mlm"

  return(res)
}

#' @export

summary.g_mlm <- function(object, digits = 3, ...) {
  varcomp <- with(object, cbind(est = c(unlist(theta), "total variance" = r_theta),
                                se = c(unlist(SE_theta), r_theta * sqrt(2 / nu))))
  betas <- with(object, cbind(est = c(coefficients$fixed, "treatment effect at a specified time" = p_beta),
                              se = c(sqrt(diag(varFix)), kappa * sqrt(r_theta))))
  ES <- with(object, cbind(est = c("unadjusted effect size" = delta_AB, "adjusted effect size" = g_AB,
                                   "degree of freedom" = nu, "constant kappa" = kappa, logLik = logLik),
                           se = c(SE_g_AB / J_nu, SE_g_AB, NA, NA, NA)))

  round(rbind(varcomp, betas, ES), digits)
}

#' @export

print.g_mlm <- function(x, digits = 3, ...) {
  ES <- with(x, cbind(est = c("unadjusted effect size" = delta_AB,
                              "adjusted effect size" = g_AB,
                              "degree of freedom" = nu),
                      se = c(SE_g_AB / J_nu, SE_g_AB, NA)))
  res <- data.frame(round(ES, digits))
  res[is.na(res)] <- ""
  res
  # print(round(ES, digits), na.print = "")
}

