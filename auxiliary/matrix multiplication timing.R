library(dplyr, warn.conflicts = FALSE, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(nlme, warn.conflicts = TRUE)
data(bdf, package = "mlmRev")

bdf_long <-
  bdf %>%
  pivot_longer(cols = c(IQ.verb, IQ.perf, aritPRET),
               names_to = "measure",
               values_to = "score") %>%
  select(schoolNR, pupilNR, sex, Minority, measure, score)

bdf_MVML <- lme(score ~ 0 + measure,
                random = ~ 1| schoolNR / pupilNR,
                corr = corSymm(form = ~ 1 | schoolNR / pupilNR),
                weights = varIdent(form = ~ 1 | measure),
                data = bdf_long)

mod <- bdf_MVML
struct <- mod$modelStruct$corStruct

res_Vinv_Vinv <- function(mod) {
  V_inv <- build_Sigma_mats(mod, invert = TRUE, sigma_scale = TRUE)
  Vinv_Vinv <- prod_blockblock(A = V_inv, B = V_inv)
  resids <- residuals(mod, level = 0) # get the rhat
  t(prod_matrixblock(A = t(resids), B = Vinv_Vinv)) # get the rhat*Vinv*dV
}

Vinv_Vinv_res <- function(mod) {
  V_inv <- build_Sigma_mats(mod, invert = TRUE, sigma_scale = TRUE)
  resids <- as.matrix(residuals(mod, level = 0)) # get the rhat
  Vinv_res <- prod_blockmatrix(A = V_inv, B = resids)
  prod_blockmatrix(A = V_inv, B = Vinv_res) # get the rhat*Vinv*dV
}


library(bench)
timings <- mark(res_Vinv_Vinv(mod), Vinv_Vinv_res(mod), iterations = 10)
