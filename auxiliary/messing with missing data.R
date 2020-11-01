library(nlme, quietly=TRUE, warn.conflicts=FALSE)

data(Orthodont)

Ortho_A <- lme(distance ~ age + Sex,
               data = Orthodont,
               random = ~ age | Subject)


test_after_deleting(Ortho_A, seed = 40)

set.seed(40)
mod <- Ortho_A

dat <- nlme::getData(mod)
y_var <- as.character(formula(mod)[[2]])
NA_y <- as.logical(rbinom(nrow(dat), size = 1, prob = 0.2))
dat[[y_var]][NA_y] <- NA
dat_complete <- dat[!NA_y,,drop=FALSE]

compare_omit_exclude_complete(mod, dat, NA_vals = NA_y)

x_vars <- attr(terms(formula(mod)), "term.labels")
x_vars <- unlist(lapply(x_vars, function(x) strsplit(x, ":")))
x_vars <- unique(x_vars)
NA_vals <- NA_y

for (x in x_vars) {
  NA_x <- as.logical(rbinom(nrow(dat), size = 1, prob = 0.2 / length(x_vars)))
  dat[[x]][NA_x] <- NA
  NA_vals <- NA_vals | NA_x
}
dat_complete <- dat[!NA_vals,,drop=FALSE]

compare_omit_exclude_complete(mod, dat, NA_vals)

mod_omit <- suppressWarnings(update(mod, data = dat, na.action = "na.omit"))
mod_exclude <- suppressWarnings(update(mod, data = dat, na.action = "na.exclude"))
mod_complete <- suppressWarnings(update(mod, data = dat_complete))

EI_omit <- Fisher_info(mod_omit, type = "expected")
EI_exclude <- Fisher_info(mod_exclude, type = "expected")
EI_complete <- Fisher_info(mod_complete, type = "expected")
AI_omit <- Fisher_info(mod_omit, type = "average")
AI_exclude <- Fisher_info(mod_exclude, type = "average")
AI_complete <- Fisher_info(mod_complete, type = "average")

testthat::expect_equal(EI_omit, EI_exclude)
testthat::expect_equal(AI_omit, AI_exclude)

mod <- mod_complete
type = "expected"

theta <- extract_varcomp(mod)
theta_names <- vapply(strsplit(names(unlist(theta)), split = "[.]"),
                      function(x) paste(unique(x), collapse = "."), character(1L))

r <- length(unlist(theta))

# Calculate derivative matrix-lists

dV_list <- build_dV_list.lme(mod)

# block-diagonal V^-1
V_inv <- build_Sigma_mats.lme(mod, invert = TRUE, sigma_scale = TRUE)

invert <- TRUE
sigma_scale <- TRUE


# lowest-level covariance structure
V_list <- build_var_cor_mats(mod, sigma_scale = sigma_scale)

# random effects covariance structure
ZDZ_list <- build_RE_mats(mod, sigma_scale = sigma_scale)

V_grps <- attr(V_list, "groups")

# Check if lowest-level covariance structure is nested within RE structure
ZDZ_grps <- attr(ZDZ_list, "groups")
group_mapping <- tapply(ZDZ_grps, V_grps, function(x) length(unique(x)))
nested <- all(group_mapping == 1L)

Sigma_list <- add_bdiag(V_list, ZDZ_list, data.frame(V_grps, ZDZ_grps))
Sigma_grps <- attr(ZDZ_list, "groups")

dat <- getData(mod)
data.frame(V = sapply(V_list, nrow), ZDZ = sapply(ZDZ_list, nrow), grps = table(dat$Subject))


R_list <- build_corr_mats(mod)
