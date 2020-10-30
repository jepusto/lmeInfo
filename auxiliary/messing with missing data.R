library(nlme, quietly=TRUE, warn.conflicts=FALSE)

data(Hartnagel, package = "carData")

Hart_AR1 <- gls(fconvict ~ tfr + partic + degrees + mconvict,
                correlation=corAR1(0.3, form = ~ year), method="REML",
                data=Hartnagel)

test_after_deleting(Hart_AR1)

mod <- Hart_AR1

dat <- nlme::getData(mod)
y_var <- as.character(formula(mod)[[2]])
NA_y <- as.logical(rbinom(nrow(dat), size = 1, prob = 0.2))
dat[[y_var]][NA_y] <- NA

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

compare_omit_exclude_complete(mod, dat, NA_vals)

mod_omit <- suppressWarnings(update(mod, data = dat, na.action = "na.omit"))
mod_exclude <- suppressWarnings(update(mod, data = dat, na.action = "na.exclude"))

EI_omit <- Fisher_info(mod_omit, type = "expected")
EI_exclude <- Fisher_info(mod_exclude, type = "expected")
AI_omit <- Fisher_info(mod_omit, type = "average")
AI_exclude <- Fisher_info(mod_exclude, type = "average")

testthat::expect_equal(EI_omit, EI_exclude)
testthat::expect_equal(AI_omit, AI_exclude)



library(nlme)
data(Orthodont)

# Orthodont

Ortho_A <- lme(distance ~ age + Sex,
               data = Orthodont,
               random = ~ age | Subject)

Ortho_B_Power <- update(Ortho_A, weights = varPower()) # fitted(.) is used by default
test_after_deleting(Ortho_B_Power)

mod <- Ortho_B_Power

# NA values in response

dat <- nlme::getData(mod)
y_var <- as.character(formula(mod)[[2]])

NA_vals <- as.logical(rbinom(nrow(dat), size = 1, prob = 0.2))
dat[[y_var]][NA_vals] <- NA

compare_omit_exclude_complete(mod, dat, NA_vals)


dat_complete <- dat[!NA_vals,,drop=FALSE]

mod_omit <- suppressWarnings(stats::update(mod, data = dat, na.action = "na.omit"))
mod_exclude <- suppressWarnings(stats::update(mod, data = dat, na.action = "na.exclude"))
mod_comp <- suppressWarnings(stats::update(mod, data = dat_complete))
mod_comp$data <- dat_complete

EI_omit <- Fisher_info(mod_omit, type = "expected")
EI_exclude <- Fisher_info(mod_exclude, type = "expected")
EI_comp <- Fisher_info(mod_comp, type = "expected")
AI_omit <- Fisher_info(mod_omit, type = "average")
AI_exclude <- Fisher_info(mod_exclude, type = "average")
AI_comp <- Fisher_info(mod_comp, type = "average")


build_var_cor_mats(mod_comp, sigma_scale = FALSE)
mod <- mod_comp
R_list <- build_corr_mats(mod)
sigma_scale <- FALSE

# NA values in predictors too

x_vars <- attr(terms(formula(mod)), "term.labels")
x_vars <- unlist(lapply(x_vars, function(x) strsplit(x, ":")))
x_vars <- unique(x_vars)
NA_all <- NA_vals

for (x in x_vars) {
  NA_vals <- as.logical(rbinom(nrow(dat), size = 1, prob = 0.2 / length(x_vars)))
  dat[[x]][NA_vals] <- NA
  NA_all <- NA_all | NA_vals
}

if (inherits(mod, "gls")) dat <<- dat
compare_omit_exclude_complete(mod, dat, NA_all)
