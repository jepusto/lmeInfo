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
  select(schoolNR, pupilNR, sex, Minority, measure, score) %>%
  arrange(schoolNR, pupilNR, measure)

bdf_MVML <- lme(score ~ 0 + measure,
                random = ~ 1| schoolNR / pupilNR,
                corr = corSymm(form = ~ 1 | schoolNR / pupilNR),
                weights = varIdent(form = ~ 1 | measure),
                data = bdf_long)

# mod <- bdf_MVML
# struct <- mod$modelStruct$corStruct

# introduce random missing to bdf_long
set.seed(20200311)
bdf_long$school_id <- bdf_long$schoolNR
levels(bdf_long$school_id) <- LETTERS[1:nlevels(bdf_long$school_id)]

bdf_long_wm <-
  bdf_long %>%
  mutate(
    row_index = rbinom(n = n(), size = 1, prob = 0.9),
    score = if_else(row_index == 1, score, as.numeric(NA)),
    measure_id = as.integer(factor(measure))
  ) %>%
  filter(!is.na(score)) %>%
  select(-row_index) %>%
  arrange(schoolNR, pupilNR, measure_id)

bdf_long_shuff <-
  sample_frac(bdf_long_wm, 1) %>%
  mutate(row = row_number())

bdf_wm <- lme(score ~ 0 + measure,
              random = ~ 1| schoolNR / pupilNR,
              corr = corSymm(form = ~ measure_id | schoolNR / pupilNR),
              weights = varIdent(form = ~ 1 | measure_id),
              data = bdf_long_wm,
              control=lmeControl(msMaxIter = 100, apVar = FALSE, returnObject = TRUE))

bdf_wm_shuff <- lme(score ~ 0 + measure,
                    random = ~ 1| schoolNR / pupilNR,
                    corr = corSymm(form = ~ measure_id | schoolNR / pupilNR),
                    weights = varIdent(form = ~ 1 | measure_id),
                    data = bdf_long_shuff,
                    control=lmeControl(msMaxIter = 100, apVar = FALSE, returnObject = TRUE))

bdf_wm_id <- lme(score ~ 0 + measure,
                  random = ~ 1| school_id / pupilNR,
                  corr = corSymm(form = ~ measure_id | school_id / pupilNR),
                  weights = varIdent(form = ~ 1 | measure_id),
                  data = bdf_long_shuff,
                  control=lmeControl(msMaxIter = 100, apVar = FALSE, returnObject = TRUE))


test_that("targetVariance() works with multivariate models.", {
  test_Sigma_mats(bdf_MVML, bdf_long$schoolNR)
  test_Sigma_mats(bdf_wm, bdf_long_wm$schoolNR)
  test_Sigma_mats(bdf_wm_id, bdf_long_shuff$school_id)
  test_Sigma_mats(bdf_wm_shuff, bdf_long_shuff$school_id)
})

test_that("Derivative matrices are of correct dimension with multivariate models.", {
  test_deriv_dims(bdf_MVML)
  test_deriv_dims(bdf_wm)
  test_deriv_dims(bdf_wm_id)
  test_deriv_dims(bdf_wm_shuff)
})

test_that("Information matrices work with FIML too.", {
  test_with_FIML(bdf_MVML)
  test_with_FIML(bdf_wm)
  test_with_FIML(bdf_wm_id)
  test_with_FIML(bdf_wm_shuff)
})

test_that("Results do not depend on order of data.", {
  test_after_shuffling(bdf_MVML, seed = 26)
  test_after_shuffling(bdf_wm)
  test_after_shuffling(bdf_wm_id, seed = 22)
  test_after_shuffling(bdf_wm_shuff, seed = 21)
})

test_that("New REML calculations work.", {
  check_REML2(bdf_MVML)
  check_REML2(bdf_wm)
  check_REML2(bdf_wm_id)
  check_REML2(bdf_wm_shuff)
})

mod <- bdf_wm_shuff
struct <- mod$modelStruct$corStruct
sigma_scale <- TRUE

sigma_sq <- if (sigma_scale) mod$sigma^2 else 1
R_list <- build_corr_mats(mod)
all_groups <- (mod$groups)
sort_order <- order(do.call(order, all_groups))
all_groups[sort_order,]
sd_vec <- sqrt(sigma_sq) / nlme::varWeights(mod$modelStruct$varStruct)[sort_order]
sd_list <- split(sd_vec, attr(R_list, "groups"))
V_list <- Map(function(R, s) tcrossprod(s) * R, R = R_list, s = sd_list)

V_list <- build_var_cor_mats(mod, sigma_scale = sigma_scale)
ZDZ_list <- build_RE_mats(mod, sigma_scale = sigma_scale)
V_grps <- attr(V_list, "groups")
ZDZ_grps <- attr(ZDZ_list, "groups")
group_mapping <- tapply(ZDZ_grps, V_grps, function(x) length(unique(x)))
nested <- all(group_mapping == 1L)

sigma_sq_id <- if (sigma_scale) bdf_wm_id$sigma^2 else 1
R_list_id <- build_corr_mats(bdf_wm_id)

all_groups_id <- rev(bdf_wm_id$groups)
sort_order_id <- order(do.call(order, all_groups_id))
sd_vec_id <- sqrt(sigma_sq_id) / nlme::varWeights(bdf_wm_id$modelStruct$varStruct)[sort_order_id]
sd_list_id <- split(sd_vec_id, attr(R_list_id, "groups"))
V_list_id <- Map(function(R, s) tcrossprod(s) * R, R = R_list_id, s = sd_list_id)

V_list_id <- build_var_cor_mats(bdf_wm_id, sigma_scale = sigma_scale)
ZDZ_list_id <- build_RE_mats(bdf_wm_id, sigma_scale = sigma_scale)
V_grps_id <- attr(V_list_id, "groups")
ZDZ_grps_id <- attr(ZDZ_list_id, "groups")
group_mapping_id <- tapply(ZDZ_grps_id, V_grps_id, function(x) length(unique(x)))
nested_id <- all(group_mapping_id == 1L)

expect_equal(unblock(R_list), unblock(R_list_id))
expect_equal(unblock(V_list), unblock(V_list_id))
expect_equal(sd_vec, sd_vec_id)

Sigma_list <- build_Sigma_mats(mod, sigma_scale = TRUE)
V_list <- build_var_cor_mats(mod, sigma_scale = TRUE)
ZDZ_list <- build_RE_mats(mod, sigma_scale = TRUE)
theta <- extract_varcomp(mod)

Sigma <- unblock(Sigma_list, attr(Sigma_list, "groups"))
V_full <- unblock(V_list, attr(V_list, "groups"))
ZDZ_full <- unblock(ZDZ_list, attr(ZDZ_list, "groups"))
expect_equal(Sigma, V_full + ZDZ_full)


Tau_params <- dV_dreStruct(mod)                           # random effects structure(s)
cor_params <- dV_dcorStruct(mod)                          # correlation structure
var_params <- dV_dvarStruct(mod)                          # variance structure
sigma_sq <- dV_dsigmasq(mod)                              # sigma_sq

# Create a list of derivative matrices
dV_list <- c(unlist(Tau_params, recursive = FALSE), cor_params, var_params, sigma_sq)

# block-diagonal V^-1
V_inv <- build_Sigma_mats(mod, invert = TRUE, sigma_scale = TRUE)

# list with V^-1 dV entries
Vinv_dV <- lapply(dV_list, prod_blockblock, A = V_inv)

A <- V_inv
B <- dV_list[[3]]
