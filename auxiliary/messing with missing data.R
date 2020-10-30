library(nlme, quietly=TRUE, warn.conflicts=FALSE)

data(Orthodont)

Ortho_hom <- gls(distance ~ Subject:age + Sex,
                 data = Orthodont)


test_after_deleting(Ortho_hom, seed = 40)

set.seed(40)
mod <- Ortho_hom

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

mod <- mod_exclude
type = "expected"
