library(dplyr, warn.conflicts = FALSE, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(nlme)

skip("Not worrying about it now.")

data(bdf, package = "mlmRev")

bdf_long <-
  bdf %>%
  filter(schoolNR %in% levels(schoolNR)[1:12]) %>%
  droplevels() %>%
  pivot_longer(cols = c(IQ.verb, IQ.perf, aritPRET),
               names_to = "measure",
               values_to = "score") %>%
  select(schoolNR, pupilNR, sex, Minority, measure, score)

bdf_MVML <- lme(score ~ 0 + measure,
                random = ~ 1| schoolNR / pupilNR,
                corr = corSymm(form = ~ 1 | schoolNR / pupilNR),
                weights = varIdent(form = ~ 1 | measure),
                data = bdf_long)

# mod <- bdf_MVML
# struct <- mod$modelStruct$corStruct

test_that("targetVariance() works with multivariate models.", {
  test_Sigma_mats(bdf_MVML, bdf_long$schoolNR)
})

test_that("Derivative matrices are of correct dimension with multivariate models.", {
  test_deriv_dims(bdf_MVML)
})

test_that("Information matrices work with FIML too.", {
  test_with_FIML(bdf_MVML)
})


# introduce random missing to bdf_long

bdf_long_wm <-
  bdf_long %>%
  mutate(
    row_index = rbinom(n = n(), size = 1, prob = 0.9),
    score = if_else(row_index == 1, score, as.numeric(NA)),
    measure_id = as.integer(factor(measure))
  ) %>%
  filter(!is.na(score)) %>%
  select(-row_index)

bdf_wm <- lme(score ~ 0 + measure,
              random = ~ 1| schoolNR / pupilNR,
              corr = corSymm(form = ~ measure_id | schoolNR / pupilNR),
              weights = varIdent(form = ~ 1 | measure_id),
              data = bdf_long_wm,
              control=lmeControl(msMaxIter = 100, apVar = FALSE, returnObject = TRUE))

mod <- bdf_wm
struct <- mod$modelStruct$corStruct
sigma_scale <- TRUE

Sigma_list <- build_Sigma_mats(mod, sigma_scale = TRUE)
V_list <- build_var_cor_mats(mod, sigma_scale = TRUE)
ZDZ_list <- build_RE_mats(mod, sigma_scale = TRUE)
theta <- extract_varcomp(mod)

Sigma <- unblock(Sigma_list, attr(Sigma_list, "groups"))
V_full <- unblock(V_list, attr(V_list, "groups"))
ZDZ_full <- unblock(ZDZ_list, attr(ZDZ_list, "groups"))
expect_equal(Sigma, V_full + ZDZ_full)

bdf_long_wm %>%
  filter(schoolNR == "2") %>%
  filter(pupilNR %in% 27001:27004) %>%
  select(pupilNR, measure, measure_id)

Sigma_list[[2]][1:10,1:10]
ZDZ_list[[2]][1:10,1:10]
theta$Tau$schoolNR
sum(unlist(theta$Tau))

V_sub <- paste0("2/", 27001:27004)
V_list[V_sub]
Sigma_sub <- lapply(V_list[V_sub], function(x) x + sum(unlist(theta$Tau)))
Sigma_sub

test_that("targetVariance() works with multivariate models.", {
  test_Sigma_mats(bdf_wm, bdf_long_wm$schoolNR)
})

test_that("Derivative matrices are of correct dimension with multivariate models.", {
  test_deriv_dims(bdf_wm)
})

test_that("Information matrices work with FIML too.", {
  test_with_FIML(bdf_wm)
})
