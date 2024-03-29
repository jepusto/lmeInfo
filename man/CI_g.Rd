% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confidence-interval-functions.R
\name{CI_g}
\alias{CI_g}
\title{Calculates a confidence interval for a standardized mean difference effect size}
\usage{
CI_g(g, cover = 0.95, bound = 35, symmetric = TRUE)
}
\arguments{
\item{g}{an estimated effect size object of class \code{g_mlm}.}

\item{cover}{confidence level.}

\item{bound}{numerical tolerance for non-centrality parameter in
\code{\link[stats]{qt}}.}

\item{symmetric}{If \code{TRUE} (the default), use a symmetric confidence
interval. If \code{FALSE}, use a non-central t approximation to obtain an
asymmetric confidence interval.}
}
\value{
A vector of lower and upper confidence bounds.
}
\description{
Calculates a confidence interval for a \code{g_mlm} object,
  using either a central t distribution (for a symmetric interval) or a
  non-central t distribution (for an asymmetric interval).
}
\examples{

library(nlme)
data(Bryant2016, package = "lmeInfo")
Bryant2016_RML1 <- lme(fixed = outcome ~ treatment,
                       random = ~ 1 | school/case,
                       correlation = corAR1(0, ~ session | school/case),
                       data = Bryant2016)
Bryant2016_g1 <- g_mlm(Bryant2016_RML1, p_const = c(0,1), r_const = c(1,1,0,1),
                       infotype = "expected")
CI_g(Bryant2016_g1, symmetric = TRUE)
CI_g(Bryant2016_g1, symmetric = FALSE)

}
